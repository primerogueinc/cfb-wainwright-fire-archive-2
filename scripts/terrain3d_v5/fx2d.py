"""Flames and smoke for the 3D fire animation, composited in 2D.

The Blender render darkens burned ground and gives the front a faint glow, but
at 3-6 km a slow morning fire reads as a dark patch, not a fire. This layer
draws what is burning at each moment from the same modelled arrival times:

  - flames: every 20 m cell whose modelled arrival is within the last few
    minutes (40 min before the wind shift, when the fire crept at ~0.4 km/h;
    4 min in the run), drawn as additive glow that flickers;
  - smoke: puffs rising from those cells and drifting downwind (from the
    southwest before 3:16 p.m., the northwest after), ILLUSTRATIVE: v4's smoke
    model is not used here, and puff speed is in screen time, not the clock.

Projection is shared with troops.Scene, so flames sit on the rendered terrain.
"""
import math
import numpy as np
from PIL import Image, ImageFilter, ImageDraw, ImageChops

T_SHIFT = 15 * 60 + 27.9  # v5: the front reaches the fire at 15:27.9 (moving line fitted to 9 stations/sensors)


class FireFX:
    def __init__(self, scene, arrival_min, step=2, seed=5, flame_len=None, ros=None):
        self.S = scene
        a = np.asarray(arrival_min, np.float64)[::step, ::step]
        fl = None if flame_len is None else np.asarray(flame_len, np.float64)[::step, ::step]
        rr, cc = np.nonzero(np.isfinite(a))
        RES = scene.RES * step; H = a.shape[0]
        self.emit_scale = (step / 2.0) ** 2   # smoke puffs per area unchanged whatever the cell size
        self.cells = np.column_stack([(cc + 0.5) * RES, (H - 1 - rr + 0.5) * RES, a[rr, cc]])
        self.L = fl[rr, cc] if fl is not None else np.full(len(rr), 2.0)   # flame length, m (Byram); nan = no flames
        # how long a 20 m cell burns: the edge crosses it at the local spread rate (plus about half a minute of burn-out)
        rs = np.asarray(ros, np.float64)[::step, ::step][rr, cc] if ros is not None else np.full(len(rr), 20.0)
        self.win = np.clip(RES / np.maximum(np.nan_to_num(rs, nan=20.0), 0.5) + 0.5, 0.6, 20.0)
        post = self.cells[:, 2] >= T_SHIFT
        self.win[post] = np.minimum(self.win[post], 0.6)   # flaming zone: grass burns out in about 30 s behind the running edge
        rng = np.random.default_rng(seed)
        self.phase = rng.uniform(0, 2 * math.pi, len(self.cells))
        self.jit = rng.uniform(-0.5, 0.5, (len(self.cells), 2)) * RES
        self.z = np.array([scene.gz(x, y) for x, y, _ in self.cells])
        self.emit = rng.random(len(self.cells))  # which cells give off smoke puffs

    def _project_many(self, f, X, Y, Z):
        loc, fwd, right, up = self.S.camera(f)
        V = np.column_stack([X, Y, Z]) - loc
        z = V @ fwd
        fpx = 30.0 / 36.0 * self.S.W
        return self.S.W / 2 + fpx * (V @ right) / z, self.S.H / 2 - fpx * (V @ up) / z, z

    def draw(self, im, f, g, t):
        """Composite flames and smoke for camera frame f, global frame g, clock t (min)."""
        W, H = im.size; sx, sy = W / self.S.W, H / self.S.H; Sc = W / 960
        dt = t - self.cells[:, 2]
        win = self.win
        act = (dt >= 0) & (dt <= win) & np.isfinite(self.L)
        if not act.any():
            return im
        idx = np.nonzero(act)[0]
        inten = np.exp(-dt[idx] / (win[idx] * 0.6))
        inten *= 0.75 + 0.25 * np.sin(g * 1.7 + self.phase[idx])  # flicker
        X = self.cells[idx, 0] + self.jit[idx, 0]; Y = self.cells[idx, 1] + self.jit[idx, 1]
        Z = self.z[idx] + 2
        px, py, depth = self._project_many(f, X, Y, Z)
        ok = (depth > 0) & (px > -50) & (px < self.S.W + 50) & (py > -50) & (py < self.S.H + 50)

        # --- smoke (drawn first, under the flames) -------------------------
        wind_to = math.radians((235 if t < T_SHIFT else 306) + 180)
        wx, wy = math.sin(wind_to), math.cos(wind_to)
        smoke = Image.new('RGBA', (W, H), (0, 0, 0, 0)); sd = ImageDraw.Draw(smoke)
        AGE = 18  # frames a puff lives
        drift = 28.0 if t < T_SHIFT else 60.0  # metres per frame downwind (screen time)
        rise = 20.0                            # scene units per frame (height is x3)
        em = idx[(self.emit[idx] < (0.10 if t < T_SHIFT else 0.012) * self.emit_scale)]
        for j in em:
            x0, y0, z0 = self.cells[j, 0], self.cells[j, 1], self.z[j]
            for k in range(0, AGE, 3):
                age = (g + int(self.phase[j] * 3)) % 3 + k
                xs = x0 + wx * drift * age; ys = y0 + wy * drift * age; zs = z0 + 10 + rise * age
                qx, qy, qd = self.S.project(f, (xs, ys, zs))
                if qd <= 0:
                    continue
                r = min((3 + 0.9 * age) * Sc * 5200 / qd, 12 * Sc)
                a = int(115 * (1 - age / AGE))
                X_, Y_ = qx * sx, qy * sy
                sd.ellipse([X_ - r, Y_ - r, X_ + r, Y_ + r], fill=(205, 200, 195, a))
        smoke = smoke.filter(ImageFilter.GaussianBlur(6 * Sc))
        im = im.convert('RGBA'); im.alpha_composite(smoke)

        # --- flames: metric tongues (Byram flame length at the modelled spread rate) + bloom -----------
        glow = np.zeros((H, W, 3), np.float32); core = np.zeros((H, W, 3), np.float32)
        loc, fwd, right, up = self.S.camera(f)
        fpx = 30.0 / 36.0 * self.S.W
        sel = np.nonzero(ok)[0]
        for q in sel:
            j = idx[q]; v = float(inten[q]); d = float(depth[q])
            L = self.L[j] * (0.75 + 0.25 * math.sin(g * 2.3 + self.phase[j] * 5)) * (0.5 + 0.5 * v) * (0.8 + 0.2 * math.cos(self.phase[j] * 7))
            x, y = px[q] * sx, py[q] * sy
            hpx = L * fpx / d * sy                      # flame height on screen, true scale
            wpx = max(1.0 * Sc, 0.8 * L * fpx / d * sx)
            # bloom: the light of the burning edge, visible even where the flames themselves are sub-pixel
            r = max(2.2 * Sc, 1.6 * hpx) * (0.8 + 0.4 * v)
            x0, x1, y0, y1 = int(x - 2 * r), int(x + 2 * r) + 1, int(y - 2 * r - hpx), int(y + 2 * r) + 1
            if x1 <= 0 or y1 <= 0 or x0 >= W or y0 >= H:
                continue
            xs = np.arange(max(x0, 0), min(x1, W))[None, :]; ys = np.arange(max(y0, 0), min(y1, H))[:, None]
            dd = ((xs - x) ** 2 + (ys - (y - 0.4 * hpx)) ** 2) / (r * r)
            glow[ys, xs] += (np.exp(-dd)[..., None] * np.array([1.0, 0.42, 0.08], np.float32) * (0.55 * v))
            # tongue: narrows to the tip, white-yellow at the base to deep orange at the tip
            if hpx >= 1.0:
                ty = np.arange(max(int(y - hpx), 0), min(int(y) + 1, H))
                if len(ty):
                    qn = np.clip((y - ty) / hpx, 0, 1)[:, None]
                    half = (wpx / 2) * (1 - qn) ** 0.7 + 0.35
                    tx = np.arange(max(int(x - wpx), 0), min(int(x + wpx) + 2, W))[None, :]
                    m = np.clip(half + 0.5 - np.abs(tx - x), 0, 1)
                    col = np.stack([np.ones_like(qn), 0.92 - 0.6 * qn, 0.55 - 0.5 * qn], -1)[:, 0, :]
                    core[ty[:, None], tx] = np.maximum(core[ty[:, None], tx], m[..., None] * col[:, None, :] * (0.75 + 0.25 * v))
        import cv2
        glow = cv2.GaussianBlur(glow, (0, 0), 1.5 * Sc) + cv2.GaussianBlur(glow, (0, 0), 6 * Sc) * 0.6
        core = cv2.GaussianBlur(core, (0, 0), 0.7 * Sc) + cv2.GaussianBlur(core, (0, 0), 2.5 * Sc) * 0.9
        b_ = np.asarray(im.convert('RGB')).astype(np.float32) / 255
        out = 1 - (1 - b_) * (1 - np.clip(glow, 0, 1))        # screen blend: vivid without clipping to white
        out = np.maximum(out, np.clip(core, 0, 1)) * 0 + (out + np.clip(core, 0, 1) * (1 - out))
        base = Image.fromarray((np.clip(out, 0, 1) * 255).astype(np.uint8))
        return base
