"""Composite the full-timeline 3D reconstruction, Rubicon Ridge (working model v5) version.

Copied from part 1 (scripts/terrain3d/compose_story.py, the v4 timeline) and changed for v5:
  * fire: the weighted-median arrival field of v5's 1,440 weighted runs (layers_v5.py);
  * clock: the gust reaches CZWN (14.7 km west-northwest) at 3:16 p.m.; the fitted front (toward 114.5 deg, 77 km/h) is drawn
    crossing the scene and reaches the fire at 15:27.9, when the wind at the fire turns;
  * places: 2024 Range Standing Orders grids (Range 2, Range 11; the road head is Range 15's entrance);
  * escape scenarios (scenarios_v5.py): leaving the official Range 2 grid at the front, 5 and 10 minutes later.

Inputs: LAYERDIR (v5 layers, story.json, scenarios_v5.json), BASEDIR (the v4 story renders, reused), SPRITEDIR, SWOB csv.
Run: python3 compose_story.py LAYERDIR BASEDIR SPRITEDIR SWOB_CSV OUT.mp4
"""
import sys, os, json, math, subprocess
import numpy as np, cv2
from PIL import Image, ImageDraw, ImageFont, ImageFilter
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import troops
from troops import Scene
troops.EXAG = 1.0  # the story base renders are true scale (scene.py 'base' mode)
from fx2d import FireFX

Ld, Bd, Sd, SWOB, OUT = sys.argv[1:6]
FPS = 12
LAB = os.path.join(os.path.dirname(os.path.abspath(OUT)), os.path.basename(OUT)[:-4] + '_lab'); os.makedirs(LAB, exist_ok=True)
F = '/usr/share/fonts/opentype/inter/'
fb = lambda s: ImageFont.truetype(F + 'Inter-Bold.otf', int(s))
fm = lambda s: ImageFont.truetype(F + 'Inter-Medium.otf', int(s))
freg = lambda s: ImageFont.truetype(F + 'Inter-Regular.otf', int(s))
INK = (14, 16, 20); WHITE = (255, 255, 255); MUTED = (200, 205, 212)
AMBER = (255, 196, 64); BLUE = (110, 170, 255); RED = (255, 90, 60); GREEN = (120, 220, 140)
hm = lambda h, m=0: h * 60 + m
T_SHIFT = hm(15, 27.9)  # v5: the front reaches the fire (CZWN's gust was at 3:16)
T_CZWN = hm(15, 16)
T_FIRST = hm(12, 17)  # first detection (Landsat 8). Nothing is drawn burning before it.

P = json.load(open(os.path.join(Ld, 'points.json')))
G = P['grid']; RES = G['res']; Wm, Hm = G['w'] * RES, G['h'] * RES
arr = np.load(os.path.join(Ld, 'arrival_min.npy')).astype(np.float64)
unk = np.load(os.path.join(Ld, 'unknown.npy'))
sev = np.load(os.path.join(Ld, 'severity.npy'))
scar5 = (np.asarray(Image.open(os.path.join(Ld, 'scar.png'))) > 127).astype(np.float32)
story = json.load(open(os.path.join(Ld, os.environ.get('STORY', 'story.json'))))['frames']
TR = Scene(Ld, Sd)
FL = np.load(os.path.join(Ld, 'flame_len.npy')) if os.path.exists(os.path.join(Ld, 'flame_len.npy')) else None
RO = np.load(os.path.join(Ld, 'ros_mpm.npy')) if os.path.exists(os.path.join(Ld, 'ros_mpm.npy')) else None
FX = FireFX(TR, arr, flame_len=FL, ros=RO, step=1)
OUT_W, OUT_H = (int(v) for v in os.environ.get('OUT_RES', '1920x1080').split('x'))

# ---------------------------------------------------------------- wind (CZWN hourly)
import pandas as pd
sw = pd.read_csv(SWOB)
sw = sw[sw.time_mdt.str.startswith('2026-10-17') & sw.max_wnd_spd_pst1hr_tm.notna()]
sw = sw[sw.time_mdt < '2026-10-17 16:00']  # the 16:00 report holds the gust front; used only once the front reaches the fire
WT = np.array([int(s[11:13]) * 60 for s in sw.time_mdt]) - 30  # mid-hour of the reporting period
WD = sw.avg_wnd_dir_10m_pst10mts.values.astype(float); WS = sw.avg_wnd_spd_10m_pst10mts.values.astype(float)
WG = sw.max_wnd_spd_10m_pst1hr.values.astype(float)


def wind(t):
    if t >= T_SHIFT:  # the front at the fire (~3:28); speeds as CZWN recorded them at its own passage (82 km/h gust)
        return 307.0, 47.6, 82.0
    t = min(t, WT[-1])  # hold the last pre-front report until the gust front
    u = np.interp(t, WT, np.sin(np.radians(WD)) * WS); v = np.interp(t, WT, np.cos(np.radians(WD)) * WS)
    return (math.degrees(math.atan2(u, v)) % 360, float(np.interp(t, WT, WS)), float(np.interp(t, WT, WG)))


# ---------------------------------------------------------------- sun (night grading)
def sun_elev(t_min, lat=52.77, lon=-110.92, doy=290):
    utc = t_min / 60 + 6
    g = 2 * math.pi / 365 * (doy - 1 + (utc - 12) / 24)
    decl = (0.006918 - 0.399912 * math.cos(g) + 0.070257 * math.sin(g) - 0.006758 * math.cos(2 * g)
            + 0.000907 * math.sin(2 * g) - 0.002697 * math.cos(3 * g) + 0.00148 * math.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * math.cos(g) - 0.032077 * math.sin(g) - 0.014615 * math.cos(2 * g) - 0.040849 * math.sin(2 * g))
    ha = math.radians((utc * 60 + eqt + 4 * lon) / 4 - 180); la = math.radians(lat)
    return 90 - math.degrees(math.acos(math.sin(la) * math.sin(decl) + math.cos(la) * math.cos(decl) * math.cos(ha)))


# ---------------------------------------------------------------- escape scenarios
SC = json.load(open(os.path.join(Ld, 'scenarios_v5.json')))
RUNS = []
for r, col in zip(SC['runners'], (GREEN, AMBER, RED)):
    RUNS.append(dict(r, col=col, odds=f"{round(r['odds'] * 100):.0f}%", x0=SC['start_xy_epsg3400'][0] - G['x0'], y0=SC['start_xy_epsg3400'][1] - G['y0']))
    print(f"scenario leave {int(r['dep'] // 60)}:{int(r['dep'] % 60):02d} -> {r['fate']} after {r['dist']:.0f} m (odds {r['odds']})")
R2_ARR = 15 * 60 + SC['range2_arrival_min_after_1500']

# the fitted front (front_plane.json): t(x, y) = 15:27.9 + gradient . (x - fire centre), minutes; gradient in min/km
FPJ = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../data/rubicon-ridge/model_v5/front_plane.json')))
GX, GY = FPJ['fit']['gradient_min_per_km']
from pyproj import Transformer as _T
_fx, _fy = _T.from_crs('EPSG:4326', G['crs'], always_xy=True).transform(-110.9052, 52.7645)
FX0, FY0 = _fx - G['x0'], _fy - G['y0']


def front_line(t):
    """map points (scene coords) where the fitted front is at clock t"""
    gn = math.hypot(GX, GY); ux, uy = GX / gn, GY / gn            # unit gradient (direction of travel)
    d = (t - T_SHIFT) / gn * 1000                                  # metres from the fire centre along the travel direction
    cx, cy = FX0 + d * ux, FY0 + d * uy
    s_ = np.linspace(-9000, 9000, 120)
    X = cx - s_ * uy; Y = cy + s_ * ux
    ok = (X > 0) & (X < Wm - 1) & (Y > 0) & (Y < Hm - 1)
    return X[ok], Y[ok]


# ---------------------------------------------------------------- helpers
def lin(v):
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)


CACHE = {}


def load_key(key):
    if key not in CACHE:
        if len(CACHE) > 6:
            CACHE.pop(next(iter(CACHE)))
        bimg = Image.open(os.path.join(Bd, f'base_{key}.png')).convert('RGB')
        if bimg.size != (OUT_W, OUT_H):   # renders are 960x540; composite at full HD (text, flames and fire edge drawn sharp)
            bimg = bimg.resize((OUT_W, OUT_H), Image.LANCZOS)
        base = np.asarray(bimg).astype(np.float32) / 255
        xy = cv2.imread(os.path.join(Bd, f'xy_{key}.png'), cv2.IMREAD_UNCHANGED).astype(np.float64) / 65535
        x = lin(xy[..., 2]) * Wm; y = lin(xy[..., 1]) * Hm; ground = xy[..., 3] > 0.5
        if x.shape != (OUT_H, OUT_W):
            x = cv2.resize(x, (OUT_W, OUT_H), interpolation=cv2.INTER_LINEAR); y = cv2.resize(y, (OUT_W, OUT_H), interpolation=cv2.INTER_LINEAR)
            ground = cv2.resize(ground.astype(np.uint8), (OUT_W, OUT_H), interpolation=cv2.INTER_NEAREST) > 0
        col = x / RES - 0.5; row = (G['h'] - 1) - (y / RES - 0.5)
        A = map_coordinates(np.where(np.isfinite(arr), arr, 1e5), [row, col], order=1, mode='nearest')
        A[A > 5e4] = np.nan
        SV = map_coordinates(sev, [row, col], order=1, mode='nearest')
        UK = map_coordinates(unk.astype(np.float32), [row, col], order=1, mode='nearest')
        SC = map_coordinates(scar5, [y_ := (2 * G['h'] - 1) - (y / (RES / 2) - 0.5), x / (RES / 2) - 0.5], order=1, mode='nearest')
        CACHE[key] = dict(base=base, A=np.where(ground, A, np.nan), SV=SV, UK=UK * ground, SC=SC * ground, ground=ground)
    return CACHE[key]


def paint_fire(d, t, final_fade):
    base = d['base'].copy(); A = d['A']
    if t < T_FIRST and not final_fade:   # before the first detection: no fire drawn
        return base, np.zeros(A.shape, np.float32)
    dt = t - A
    burning = np.nan_to_num(dt, nan=-1) >= 0
    tau = np.where(A < T_SHIFT, 15.0, 0.8)
    glow = np.where(burning, np.exp(-np.clip(np.nan_to_num(dt), 0, None) / tau), 0.0) ** 2
    inscar = d['SC'] > 0.5
    char = np.where(burning & inscar, np.clip(np.nan_to_num(dt) / 2.0, 0, 1), 0.0) * (0.6 + 0.4 * d['SV'])
    charcoal = np.array([0.055, 0.045, 0.04], np.float32)
    img = base * (1 - char[..., None]) + charcoal * char[..., None]
    if final_fade > 0:
        u = (d['UK'] * 0.8 * final_fade)[..., None]
        img = img * (1 - u) + np.array([0.33, 0.30, 0.27], np.float32) * u
    return img, glow * inscar


def night(img, glow, t):
    e = sun_elev(t)
    dark = float(np.clip((2.0 - e) / 8.0, 0, 1))  # full night below -6 deg, daylight above +2
    if dark > 0:
        tint = np.array([0.78, 0.84, 1.0], np.float32)  # moonlit, readable: the night is a cue, not a blackout
        img = img * (1 - 0.5 * dark) * (tint * dark + (1 - dark))
    g = glow[..., None] * np.array([1.0, 0.36, 0.06], np.float32)
    gl = cv2.GaussianBlur(g, (0, 0), 1.2) * 1.0 + cv2.GaussianBlur(g, (0, 0), 7) * (0.45 + 1.4 * dark)
    return np.clip(img + gl, 0, 1)


def clock(t):
    h, m = int(t // 60), int(t % 60)
    return f'{h - 12 if h > 12 else h}:{m:02d} ' + ('p.m.' if h >= 12 else 'a.m.')


def caption(fr, t):
    if fr['final']:
        return ('Oct. 24', 'The burn scar, as Sentinel-2 saw it a week later',
                'Grey: burned, timing unknown (not connected to the main fire in the model; includes heat VIIRS saw on Oct. 14).', WHITE)
    if fr['card'] == 'alert':
        return (clock(t), 'Environment Canada issues a wind warning for areas to the west',
                'Its 4:30 a.m. forecast had already called for west winds gusting to 80 km/h this afternoon. No fire has been detected yet.', RED)
    if t < T_FIRST:
        return (clock(t), 'No fire detected yet',
                'First detection: Landsat 8, 12:17 p.m. The AGL was in use on Range 11 (Maj.-Gen. Rutland); positions illustrative.', WHITE)
    if t < hm(12, 18):
        return (clock(t), 'Landsat 8: the first detection of the fire',
                'Yellow: 59 fire pixels around a strip about 2 km long that had already burned: burning for some time by then (inferred).', AMBER)
    if fr['card'] == 'warning':
        return (clock(t), 'Environment Canada adds Wainwright to its wind warning: gusts up to 100 km/h',
                'Wainwright added to the wind warning in force to the west since 11:42 a.m. The fire is about 1.2 km from Range 2.', RED)
    if t < hm(14, 49):
        return (clock(t), 'The fire creeps north; the wind rises', 'Spread between satellite passes is modelled (Rubicon Ridge v5).', AMBER)
    if t < T_CZWN:
        return (clock(t), 'NOAA-20 sees the fire (blue)', 'The wind is gusting to 50-55 km/h from the southwest.', BLUE)
    if fr['card'] == 'czwn':
        return (clock(t), 'An 82 km/h gust at the CFB Wainwright station (CZWN), 14.7 km west-northwest',
                'The cold front is still on its way: moving east-southeast at about 77 km/h, it reaches the fire about 12 minutes later.', AMBER)
    if t < T_SHIFT:
        return (clock(t), 'The front crosses the range (dashed: fitted to 9 stations and sensors)',
                'At the fire the wind is still from the southwest.', AMBER)
    if fr['card'] == 'front':
        return ('3:28 p.m.', 'The front reaches the fire. The wind turns northwest',
                'The fire\'s long southeast side becomes its front. Its edge is about 910 m from Range 2 and 1.1 km from Range 15 / the road head.', RED)
    if t < R2_ARR:
        return (clock(t), 'Three scenarios: leaving Range 2 at 3:28, 3:33 or 3:38',
                'Jogging in kit, best route: east, about 200 m to ground the fire never reaches (v5 median fire). Scenarios, not what happened.', RED)
    if t < hm(15, 46):
        return (clock(t), 'The fire reaches Range 2 at about 3:39 and Range 15 / the road head at about 3:42',
                'Weighted medians of v5\'s runs (10-90%: 3:34-3:51 and 3:36-3:58). Measured from the front, close to v4.', RED)
    return (clock(t), 'The fire runs on southeast', 'The southeast lobe burns until after 4 p.m. (modelled).', RED)


def arrow(d, cx, cy, bearing_screen, length, col, width):
    a = math.radians(bearing_screen); dx, dy = math.sin(a), -math.cos(a)
    x0, y0 = cx - dx * length / 2, cy - dy * length / 2; x1, y1 = cx + dx * length / 2, cy + dy * length / 2
    d.line([x0, y0, x1, y1], fill=col, width=width); hx, hy = -dy, dx
    d.polygon([(x1 + dx * 10, y1 + dy * 10), (x1 - dx * 8 + hx * 9, y1 - dy * 8 + hy * 9),
               (x1 - dx * 8 - hx * 9, y1 - dy * 8 - hy * 9)], fill=col)


def label_at(d, x, y, text, sub, S, anchor='la', col=WHITE):
    d.text((x, y), text, font=fb(14 * S), fill=col, stroke_width=3, stroke_fill=INK, anchor=anchor)
    if sub:
        d.text((x, y + 17 * S), sub, font=freg(11 * S), fill=MUTED, stroke_width=3, stroke_fill=INK, anchor=anchor)


def runner_state(r, t):
    if t < r['dep']:
        return None
    dist = min((t - r['dep']) * 8000 / 60, r['dist'])
    done = t >= r['t_end']
    b = math.radians(r['dir'])
    return r['x0'] + dist * math.sin(b), r['y0'] + dist * math.cos(b), dist, done


# ---------------------------------------------------------------- main loop
N = len(story)
ONLY = os.environ.get('ONLY')  # e.g. "40,120,300" to compose a few frames for checking
sel = set(int(v) for v in ONLY.split(',')) if ONLY else None
for fr in story:
    if sel is not None and fr['i'] not in sel:
        continue
    i, t = fr['i'], fr['t']
    TR.override = (fr['loc'], fr['look'])
    d = load_key(fr['key'])
    final_fade = 1.0 if fr['final'] else 0.0
    img, glow = paint_fire(d, t, final_fade)
    img = night(img, glow, t)
    im = Image.fromarray((img * 255).astype(np.uint8))
    W, H = im.size; S = W / 960
    e = sun_elev(t); dark = float(np.clip((2.0 - e) / 8.0, 0, 1))
    if not fr['final'] and t >= T_FIRST:
        im = FX.draw(im, i, i, t)  # flames and smoke (none before the first detection)
    # grenadiers: on the firing line (illustrative)
    if not fr['final'] and t >= hm(7, 50):
        TR.draw(im, i, i, min(t, T_SHIFT - 1), group=TR.r11, alpha=1.0)
        if t < hm(12, 0):
            TR.draw_agl(im, i, i % 48, 48)
    # C Company: moving until 3:16; dimmed while the scenarios play
    if not fr['final'] and t >= hm(7, 50):
        a_c = 1.0 if t < T_SHIFT else 0.35
        TR.draw(im, i, i, t, alpha=a_c * (1 - 0.6 * dark))
    dr = ImageDraw.Draw(im)
    # satellite detections
    for hs in P['hotspots']:
        if '2026-10-17' not in (hs['time_mdt'] or ''):
            continue
        if hs['source'] == 'LANDSAT_NRT' and hm(12, 17) <= t < T_SHIFT:
            col, rr = AMBER, 2.5
        elif hs['source'] == 'VIIRS_NOAA20_NRT' and hm(14, 49) <= t < T_SHIFT:
            col, rr = BLUE, 5
        else:
            continue
        x, y, z = TR.project(i, (hs['x'], hs['y'], TR.gz(hs['x'], hs['y'])))
        if z > 0:
            x *= W / TR.W; y *= H / TR.H; q = rr * S
            dr.rectangle([x - q, y - q, x + q, y + q], outline=col, width=max(2, int(2 * S)))
    # the fitted front crossing the scene (CZWN gust 3:16 -> the fire 3:28 -> on east-southeast)
    if T_CZWN <= t <= T_SHIFT + 4 and not fr['final']:
        FXs, FYs = front_line(t)
        if len(FXs) > 1:
            q = [TR.project(i, (x_, y_, TR.gz(x_, y_))) for x_, y_ in zip(FXs, FYs)]
            q = [(a_ * W / TR.W, b_ * H / TR.H) for a_, b_, c_ in q if c_ > 0]
            for k_ in range(0, len(q) - 1, 2):
                dr.line([q[k_], q[k_ + 1]], fill=INK, width=int(6 * S)); dr.line([q[k_], q[k_ + 1]], fill=BLUE, width=int(3 * S))
            if q:
                qq = min(q, key=lambda p_: p_[1] if 110 * S < p_[1] else 1e9)
                label_at(dr, qq[0] + 8 * S, qq[1] + 4 * S, 'cold front (fitted)', None, S, col=BLUE)
    # scenario runners
    if T_SHIFT <= t < hm(15, 46) and not fr['final']:
        for jr, r in enumerate(RUNS):
            st = runner_state(r, t)
            if st is None:
                continue
            x, y, dist, done = st
            pts = []
            for s_ in np.linspace(0, dist, 20):
                b = math.radians(r['dir']); px_, py_ = r['x0'] + s_ * math.sin(b), r['y0'] + s_ * math.cos(b)
                q = TR.project(i, (px_, py_, TR.gz(px_, py_))); pts.append((q[0] * W / TR.W, q[1] * H / TR.H))
            if len(pts) > 1:
                dr.line(pts, fill=INK, width=int(5 * S)); dr.line(pts, fill=r['col'], width=int(3 * S))
            sx, sy, sz = TR.project(i, (x, y, TR.gz(x, y)))
            sx *= W / TR.W; sy *= H / TR.H
            key = ('standing_observer' if done else 'walking_equipment') + '_' + ('right' if pts[-1][0] >= pts[0][0] else 'left')
            sp = TR.sprites[key]; k = 22 * S * min(5200 / sz, 1.6) * 0.75 / sp.height
            spr = sp.resize((max(1, int(sp.width * k)), max(1, int(sp.height * k))), Image.LANCZOS)
            alpha = 1.0 if not (done and r['fate'] == 'overtaken') else max(0.3, 1 - (t - r['t_end']) / 1.5)
            a_ = spr.split()[3].point(lambda v: int(v * alpha)); halo = a_.filter(ImageFilter.MaxFilter(3))
            dark_ = Image.new('RGBA', spr.size, (10, 10, 12, 0)); dark_.putalpha(halo)
            body = spr.copy(); body.putalpha(a_)
            pos = (int(sx - spr.width / 2), int(sy - spr.height)); im.paste(dark_, pos, dark_); im.paste(body, pos, body)
            dr = ImageDraw.Draw(im)
            tag = f"leave {clock(r['dep'])[:4]}"
            if done:
                tag += {'clear': f" · clear {clock(r['t_end'])[:4]}",
                        'close': f" · clear {clock(r['t_end'])[:4]}, fire under 3 min behind",
                        'overtaken': (" · fire reached Range 2 first" if r['dist'] == 0
                                      else f" · overtaken {clock(r['t_end'])[:4]}")}[r['fate']]
            lx, anc = (sx + 10 * S, 'la') if sx < W - 420 * S else (sx - 10 * S, 'ra')   # keep clear of the wind compass
            label_at(dr, lx, sy - (30 + 17 * (2 - jr)) * S, tag, None, S, anchor=anc, col=r['col'])
        # legend
        bx, by = 18 * S, 108 * S
        box = Image.new('RGBA', (int(360 * S), int(92 * S)), (*INK, 205)); im.paste(box, (int(bx), int(by)), box)
        dr = ImageDraw.Draw(im)
        dr.text((bx + 10 * S, by + 8 * S), 'Got clear (3-min margin), share of v5\'s weighted runs', font=fb(12 * S), fill=WHITE)
        for j, r in enumerate(RUNS):
            dr.rectangle([bx + 10 * S, by + (30 + 19 * j) * S, bx + 24 * S, by + (42 + 19 * j) * S], fill=r['col'])
            dr.text((bx + 32 * S, by + (28 + 19 * j) * S), f"leave {clock(r['dep'])}: {r['odds']}", font=fm(12 * S), fill=WHITE)
    # place labels
    lab_pts = {'Range 11': P['extra']['Range 11']}
    v = P['extra']['Range 11']; x, y, z = TR.project(i, (v['x'], v['y'], TR.gz(v['x'], v['y'])))
    x *= W / TR.W; y *= H / TR.H
    if z > 0 and x < 0 and 100 * S < y < H - 50 * S and not fr['final']:
        dr.text((10 * S, y), '\u2190 Range 11 (AGL range)', font=fb(14 * S), fill=WHITE, stroke_width=3, stroke_fill=INK, anchor='lm')
    if t >= T_SHIFT - 1:
        lab_pts.update({k: P['pov'][k] for k in ('Range 2', 'Roadhead')})
    NICE = {'Roadhead': 'Range 15 / road head'}
    for k, v in lab_pts.items():
        x, y, z = TR.project(i, (v['x'], v['y'], TR.gz(v['x'], v['y'])))
        x *= W / TR.W; y *= H / TR.H
        if z > 0 and 0 < x < W and 100 * S < y < H - 50 * S:
            dr.ellipse([x - 4 * S, y - 4 * S, x + 4 * S, y + 4 * S], fill=WHITE, outline=INK, width=2)
            anc, off = (('rb', (-8, -6)) if k == 'Range 2' else ('lt', (8, 3)))
            if anc == 'lt' and x > W - 230 * S:
                anc, off = 'rt', (-8, 3)
            dr.text((x + off[0] * S, y + off[1] * S), NICE.get(k, k) if t >= T_SHIFT - 1 else k, font=fb(15 * S), fill=WHITE, stroke_width=3, stroke_fill=INK, anchor=anc)
    # header
    hdr = Image.new('RGBA', (W, int(92 * S)), (*INK, 210)); im.paste(hdr, (0, 0), hdr); dr = ImageDraw.Draw(im)
    big, cap, sub, col = caption(fr, t)
    dr.text((18 * S, 12 * S), big, font=fb(34 * S), fill=col)
    cx = 18 * S + dr.textlength(big, font=fb(34 * S)) + 18 * S
    dr.text((cx, 16 * S), cap, font=fb(17 * S), fill=WHITE)
    dr.text((cx, 44 * S), sub, font=freg(12 * S), fill=MUTED)
    dr.text((W - 14 * S, 72 * S), 'RUBICON RIDGE · WORKING MODEL v5 · WORKING RESULTS, NOT YET PUBLISHED', font=fb(10 * S), fill=AMBER, anchor='rm')
    # wind compass (screen-relative)
    if not fr['final']:
        loc, look = np.array(fr['loc']), np.array(fr['look'])
        heading = math.degrees(math.atan2(look[0] - loc[0], look[1] - loc[1])) % 360
        wd, ws, wg = wind(t)
        bx, by = W - 70 * S, 150 * S
        dr.ellipse([bx - 46 * S, by - 46 * S, bx + 46 * S, by + 46 * S], fill=INK, outline=MUTED, width=2)
        arrow(dr, bx, by, (wd + 180 - heading) % 360, 60 * S, RED if wg >= 70 else (AMBER if wg >= 40 else WHITE), int(5 * S))
        nx, ny = bx + 34 * S * math.sin(math.radians(-heading)), by - 34 * S * math.cos(math.radians(-heading))
        dr.text((nx, ny), 'N', font=fb(13 * S), fill=MUTED, anchor='mm')
        dr.text((bx, by + 60 * S), f'Wind {ws:.0f} km/h', font=fm(13 * S), fill=WHITE, anchor='mm', stroke_width=2, stroke_fill=INK)
        dr.text((bx, by + 76 * S), f'gusts {wg:.0f}', font=fm(13 * S), fill=WHITE, anchor='mm', stroke_width=2, stroke_fill=INK)
    # footer
    ft = Image.new('RGBA', (W, int(40 * S)), (*INK, 255)); im.paste(ft, (0, H - int(40 * S)), ft); dr = ImageDraw.Draw(im)
    dr.text((18 * S, H - 33 * S), 'Clock not to scale · true scale · fire: Rubicon Ridge working model v5, weighted median of 1,440 runs (modelled) · '
            'observed: Landsat 8, VIIRS, CZWN wind, ECCC, Sentinel-2 scar', font=freg(11 * S), fill=MUTED)
    dr.text((18 * S, H - 17 * S), 'Figures representative · flames true scale (Byram, modelled spread) · smoke illustrative · ranges: 2024 Standing Orders · '
            'working results, not yet published · Prime Rogue Inc.', font=freg(11 * S), fill=MUTED)
    im.convert('RGB').save(os.path.join(LAB, f'{i:04d}.png'))
    if i % 40 == 0:
        print('frame', i, clock(t), flush=True)

if sel is not None:
    sys.exit(0)
subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(FPS), '-i', os.path.join(LAB, '%04d.png'),
                '-vf', 'minterpolate=fps=24:mi_mode=blend', '-c:v', 'libx264', '-preset', 'veryfast', '-crf', '19',
                '-pix_fmt', 'yuv420p', '-movflags', '+faststart', OUT], check=True)
print('wrote', OUT, N, 'frames')
