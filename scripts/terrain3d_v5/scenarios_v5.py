"""Escape scenarios for the v5 3D timeline: leaving the official Range 2 grid (and the road head) at the front's arrival at the
fire, 3 minutes later and 7 minutes later.

Method as v4 part 4 (part 1, data/authentic-stoic-stress-test/part4_escape/), simplified: straight routes in 72 directions,
jogging in kit at 8 km/h on flat ground (no slope adjustment), safe ground = never reached by the run (outside the scar + 150 m
confinement) or burned at least 5 minutes earlier; a route works if it reaches safe ground with the fire at least 3 minutes
behind at every point. Departures are measured from each run's own front time at the fire (15:27, 15:29 or 15:31).

Odds: share of the fit weight (1,440 weighted runs, v5 weighting) with at least one working route. The video's runners use the
weighted-median arrival field and the route that works best when leaving at the front (shortest distance to safety).

Run: python3 scenarios_v5.py GRIDS_DIR (rerun grids_*.npz + v5_median.npz) OUT.json
"""
import sys, os, glob, json, math
import numpy as np, pandas as pd
from pyproj import Transformer

GD, OUT = sys.argv[1:3]
R2 = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))
z = np.load(os.path.join(GD, 'v5_median.npz'))
minx, maxy, RES, H, W = float(z['minx']), float(z['maxy']), float(z['res']), int(z['H']), int(z['W'])
tr = Transformer.from_crs(4326, str(z['crs']), always_xy=True)
PTS = {'Range 2 (RSO)': (-110.90089, 52.772), 'Road head': (-110.89965, 52.77045)}
DS = np.arange(0, 4001, 20.0); V = 8000 / 60; BRG = np.arange(0, 360, 5)
DELAYS = tuple(int(v) for v in os.environ.get("DELAYS", "0,5,10").split(","))


def routes(a, x0, y0):
    """fire arrival (min after 15:00, nan = never) along every route: array [bearing, distance]"""
    X = x0 + DS[None] * np.sin(np.radians(BRG))[:, None]; Y = y0 + DS[None] * np.cos(np.radians(BRG))[:, None]
    c = ((X - minx) // RES).astype(int); r = ((maxy - Y) // RES).astype(int)
    inside = (r >= 0) & (r < H) & (c >= 0) & (c < W)
    return np.where(inside, a[np.clip(r, 0, H - 1), np.clip(c, 0, W - 1)], np.nan)


def evaluate(fa, dep):
    """per bearing: (works, index of first safe point, fate, end index) for departure dep (min after 15:00)"""
    t = dep + DS / V
    safe = np.isnan(fa) | ((fa <= t[None] - 5) & (np.arange(len(DS))[None] > 0))
    fire3 = (fa <= t[None] + 3) & ~safe; fire0 = (fa <= t[None]) & ~safe
    first = lambda m: np.where(m.any(1), m.argmax(1), 10 ** 9)
    i_s, i_3, i_0 = first(safe), first(fire3), first(fire0)
    return i_s, i_3, i_0


res = {}
names, grids = [], []
for f in sorted(glob.glob(os.path.join(GD, 'grids_*.npz'))):
    q = np.load(f); names += list(q['runs']); grids.append(q['grids'])
G = np.concatenate(grids).astype(np.float32); del grids
G[G >= 210 - 1e-6] = np.nan
df = pd.read_csv(os.path.join(R2, 'data/rubicon-ridge/model_v5/ensemble_runs.csv')).set_index('run').loc[names]
w = df.weight.values / df.weight.values.sum()
tfront = df.front.map({'15:27': 27.0, '15:29': 29.0, '15:31': 31.0}).values
for k, (lo, la) in PTS.items():
    x0, y0 = tr.transform(lo, la)
    ok = np.zeros((len(names), len(DELAYS)), bool)
    for n in range(len(names)):
        fa = routes(G[n], x0, y0)
        for j, d in enumerate(DELAYS):
            i_s, i_3, _ = evaluate(fa, tfront[n] + d)
            ok[n, j] = bool(((i_s < i_3) & (i_s < 10 ** 9)).any())
    res[k] = {f'front+{d}': round(float((ok[:, j] * w).sum()), 3) for j, d in enumerate(DELAYS)}
    print(k, res[k])

# runners on the median field
med = z['median']; T0 = 28.0   # front at the fire 15:27.9, rounded to the minute for the runners
x0, y0 = tr.transform(*PTS['Range 2 (RSO)'])
fa = routes(med, x0, y0)
i_s, i_3, i_0 = evaluate(fa, T0)
good = np.nonzero((i_s < i_3) & (i_s < 10 ** 9))[0]
b = int(good[np.argmin(DS[i_s[good]])]) if len(good) else int(np.argmin(i_s))
runs = []
for d in DELAYS:
    i_s, i_3, i_0 = evaluate(fa, T0 + d)
    s, f3, f0 = int(i_s[b]), int(i_3[b]), int(i_0[b])
    if s < f3: end, fate = s, 'clear'
    elif s < f0: end, fate = s, 'close'
    else: end, fate = f0, 'overtaken'
    t_end = 15 * 60 + T0 + d + DS[end] / V
    runs.append(dict(dep=15 * 60 + T0 + d, delay=d, dir=int(BRG[b]), dist=float(DS[end]), t_end=float(t_end), fate=fate,
                     odds=res['Range 2 (RSO)'][f'front+{d}']))
    print(f'median field: leave front+{d} -> {fate} after {DS[end]:.0f} m at {int(t_end // 60)}:{t_end % 60:04.1f} (bearing {BRG[b]})')
r2 = med[int((maxy - y0) // RES), int((x0 - minx) // RES)]
json.dump(dict(method=__doc__.split('Run:')[0].strip(), odds=res, runners=runs, range2_arrival_min_after_1500=float(r2),
               start_xy_epsg3400=[x0, y0]), open(OUT, 'w'), indent=1)
