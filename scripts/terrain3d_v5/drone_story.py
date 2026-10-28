"""Storyboard for a hypothetical drone that follows the fire (Rubicon Ridge v5), for compose_story.py.

Nothing is drawn burning before 12:17 p.m., the first satellite detection. The drone holds a wide view from the
11:42 a.m. wind warning until the front reaches the fire (3:28), then descends and follows the fire's edge: first the
edge nearest Range 2 (it reaches Range 2 at about 3:39 on v5's weighted median), then the head of the run southeast,
and pulls back up at the end. Its path is derived from the v5 weighted-median arrival field (arrival_min.npy).

Moving shots are rendered at 6 camera positions per second (every second frame); compose_story.py draws 12 frames per
second and the encoder blends between them.

Run: python3 drone_story.py LAYERDIR   (writes LAYERDIR/story_drone.json)
"""
import sys, os, json, math
import numpy as np
from scipy.ndimage import gaussian_filter1d

FPS = 12
hm = lambda h, m=0: h * 60 + m
Ld = sys.argv[1]
P = json.load(open(os.path.join(Ld, 'points.json')))
h = np.load(os.path.join(Ld, 'height.npy')).astype(np.float64); Z0 = float(np.percentile(h, 1))
arr = np.load(os.path.join(Ld, 'arrival_min.npy')).astype(np.float64)
RES = P['grid']['res']; H, W = arr.shape
rr, cc = np.nonzero(np.isfinite(arr)); CX = (cc + 0.5) * RES; CY = (H - 1 - rr + 0.5) * RES; CT = arr[rr, cc]
T_FRONT = hm(15, 27.9)


def gz(x, y):
    c = int(np.clip(x / RES, 0, W - 1)); r = int(np.clip(H - 1 - y / RES, 0, H - 1)); return h[r, c] - Z0


r2 = P['pov']['Range 2']


def target(t):
    m = (CT > t - 1.5) & (CT <= t)
    if not m.any():
        m = (CT <= t) & (CT > t - 6)
    x, y = CX[m], CY[m]
    if t < hm(15, 41):                      # the edge nearest Range 2
        k = np.argmin(np.hypot(x - r2['x'], y - r2['y']))
    else:                                    # the head of the run (farthest southeast)
        k = np.argmax(x * math.sin(math.radians(135)) + y * math.cos(math.radians(135)))
    return x[k], y[k]


# wide pose over the whole pre-front fire (from the south-southwest)
pre = CT <= T_FRONT
cx, cy = CX[pre].mean(), CY[pre].mean()
ang = math.radians(205)
WIDE = (np.array([cx + 3800 * math.sin(ang), cy + 3800 * math.cos(ang), gz(cx, cy) + 2200]), np.array([cx + 300, cy, gz(cx, cy)]))

# moving segment: clock knots and camera distance / height (zoom) knots
MOVE = [  # (screen seconds, clock)
    (0.0, T_FRONT), (6.0, hm(15, 34)), (14.0, hm(15, 39.5)), (19.0, hm(15, 43)), (27.0, hm(15, 58)), (33.0, hm(16, 15)), (36.0, hm(16, 20))]
ZOOM = [  # (screen seconds, back distance m, height m)
    (0.0, 3800, 2200), (5.0, 1500, 800), (11.0, 650, 330), (17.0, 600, 300), (22.0, 1100, 550), (31.0, 1300, 650), (36.0, 3600, 2100)]
DUR = MOVE[-1][0]
n_key = int(DUR * FPS / 2) + 1
ks = np.arange(n_key) * 2 / FPS
tk = np.interp(ks, *zip(*MOVE))
tg = np.array([target(t) for t in tk])
tg = np.column_stack([gaussian_filter1d(tg[:, 0], 3, mode='nearest'), gaussian_filter1d(tg[:, 1], 3, mode='nearest')])
# blend the first keys from the wide view's look point into the tracked edge
w0 = np.clip(ks / 5.0, 0, 1); w0 = w0 * w0 * (3 - 2 * w0)
w1 = np.clip((ks - 31.0) / 5.0, 0, 1); w1 = w1 * w1 * (3 - 2 * w1)
look_x = WIDE[1][0] * (1 - w0) + tg[:, 0] * w0; look_y = WIDE[1][1] * (1 - w0) + tg[:, 1] * w0
look_x = look_x * (1 - w1) + WIDE[1][0] * w1; look_y = look_y * (1 - w1) + WIDE[1][1] * w1
dist = np.interp(ks, [z[0] for z in ZOOM], [z[1] for z in ZOOM]); hgt = np.interp(ks, [z[0] for z in ZOOM], [z[2] for z in ZOOM])
bear = math.radians(205)
keys = []
for k in range(n_key):
    lz = gz(look_x[k], look_y[k])
    loc = [look_x[k] + dist[k] * math.sin(bear), look_y[k] + dist[k] * math.cos(bear), 0.0]
    loc[2] = max(lz + hgt[k], gz(loc[0], loc[1]) + 80)
    keys.append((loc, [look_x[k], look_y[k], lz]))
# first key equals the wide pose
keys[0] = (WIDE[0].tolist(), WIDE[1].tolist())

frames = []


def add(dur, t0, t1, key, loc, look, card=None, final=False):
    n = int(round(dur * FPS))
    for i in range(n):
        frames.append(dict(i=len(frames), seg=0, q=i / n, t=t0 + (t1 - t0) * i / n, key=key, loc=list(loc), look=list(look),
                           card=card, final=final))


wl, wk = WIDE[0].tolist(), WIDE[1].tolist()
add(4.0, hm(11, 42), hm(11, 42), 'DW', wl, wk, card='alert')
add(3.0, hm(11, 42), hm(12, 17), 'DW', wl, wk)
add(3.0, hm(12, 17), hm(12, 17), 'DW', wl, wk)
add(2.5, hm(12, 17), hm(14, 14), 'DW', wl, wk)
add(3.5, hm(14, 14), hm(14, 18), 'DW', wl, wk, card='warning')
add(2.0, hm(14, 18), hm(15, 16), 'DW', wl, wk)
add(3.5, hm(15, 16), hm(15, 16.2), 'DW', wl, wk, card='czwn')
add(2.5, hm(15, 16.2), T_FRONT, 'DW', wl, wk)
add(3.0, T_FRONT, T_FRONT + 0.1, 'DW', wl, wk, card='front')
for i in range(int(DUR * FPS)):
    s = i / FPS; k = min(i // 2, n_key - 1)
    loc, look = keys[k]
    frames.append(dict(i=len(frames), seg=1, q=s / DUR, t=float(np.interp(s, *zip(*MOVE))), key=f'D_{k:03d}', loc=loc, look=look,
                       card=None, final=False))
add(4.0, hm(16, 20), hm(16, 20), f'D_{n_key - 1:03d}', *keys[-1], final=True)
json.dump({'frames': frames}, open(os.path.join(Ld, 'story_drone.json'), 'w'))
print(len(frames), 'frames', round(len(frames) / FPS, 1), 's;', len({f['key'] for f in frames}), 'renders')
