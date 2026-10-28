"""Step 2 after v5_rerun_grids.py. Check the re-run against the published ensemble_runs.csv, then build the weighted-median arrival grid (+ p10/p90, burn
probability) and the best run's grid. Outputs in this folder only."""
import os, sys, json, glob
import numpy as np, pandas as pd
from pyproj import Transformer
O = os.path.dirname(os.path.abspath(__file__))
R2 = "/home/claude/cfb-wainwright-fire-archive-2"
meta = json.load(open(O + "/grid_meta.json"))
names, grids = [], []
for f in sorted(glob.glob(O + "/grids_*.npz")):
    z = np.load(f); names += list(z["runs"]); grids.append(z["grids"])
G = np.concatenate(grids).astype("float32"); del grids; G[G >= 210.0 - 1e-9] = np.inf   # past t_end (210 min) = not reached
G[~np.isfinite(G)] = np.inf
df = pd.read_csv(R2 + "/data/rubicon-ridge/model_v5/ensemble_runs.csv").set_index("run").loc[names]
w = df.weight.values; w = w / w.sum()
tr = Transformer.from_crs(4326, meta["crs"], always_xy=True)
def cell(lon, lat):
    x, y = tr.transform(lon, lat); return int((meta["maxy"] - y) // meta["res"]), int((x - meta["minx"]) // meta["res"])
pts = {"Range 2 (RSO)": (-110.90089, 52.772), "Range 15 (RSO)": (-110.89946, 52.77034), "Road head": (-110.89965, 52.77045),
       "Range 2 point (OSM, v4)": (-110.9063, 52.7728)}
# 1) reproduction check against the published per-run arrivals (stored at 0.01 min)
for k in ["Range 2 (RSO)", "Range 15 (RSO)", "Road head", "Range 2 point (OSM, v4)"]:
    i, j = cell(*pts[k]); mine = G[:, i, j]; pub = df[f"arr_scar|{k}"].values
    ok = np.isfinite(pub); d = np.abs(np.where(np.isfinite(mine), mine, 1e9)[ok] - pub[ok])
    print(f"{k}: max |diff| {d.max():.3f} min over {ok.sum()} runs")
# 2) per-cell weighted quantiles (inf = unreached)
def wq(a, q):
    o = np.argsort(a, axis=0); s = np.take_along_axis(a, o, 0); ws = w[o]
    c = np.cumsum(ws, axis=0) - 0.5 * ws
    idx = (c < q).sum(0); idx = np.clip(idx, 0, a.shape[0] - 1)
    return np.take_along_axis(s, idx[None], 0)[0]
med = np.empty(G.shape[1:]); p10 = np.empty_like(med); p90 = np.empty_like(med)
for r0 in range(0, G.shape[1], 20):
    blk = G[:, r0:r0 + 20]; med[r0:r0 + 20] = wq(blk, 0.5); p10[r0:r0 + 20] = wq(blk, 0.1); p90[r0:r0 + 20] = wq(blk, 0.9)
pburn = (np.isfinite(G) * w[:, None, None]).sum(0)
best = names[int(np.argmax(df.weight.values))]
gbest = G[names.index(best)]
for k, v in pts.items():
    i, j = cell(*v)
    f = lambda t: "-" if not np.isfinite(t) else f"{int(15*60+t)//60}:{int(round(15*60+t))%60:02d}"
    print(f"{k:26s} median {f(med[i,j])}  p10 {f(p10[i,j])}  p90 {f(p90[i,j])}  best run {f(gbest[i,j])}  P(burn) {pburn[i,j]:.2f}")
print("best run:", best)
np.savez_compressed(O + "/v5_median.npz", median=np.where(np.isfinite(med), med, np.nan), p10=np.where(np.isfinite(p10), p10, np.nan),
                    p90=np.where(np.isfinite(p90), p90, np.nan), pburn=pburn, best=np.where(np.isfinite(gbest), gbest, np.nan),
                    best_run=best, **{k: v for k, v in meta.items() if k != "crs"}, crs=meta["crs"])
