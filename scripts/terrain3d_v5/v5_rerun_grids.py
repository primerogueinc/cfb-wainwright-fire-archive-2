"""Re-run the 1,440 weighted v5 scar-limited solves (published confinement: revised scar + 150 m), as v5_ensemble.py
lines 101-125, and keep full arrival grids (not stored in the repository: about 0.7 GB).
Writes only next to this script's copy in a scratch folder; run it from there. Usage: v5_rerun_grids.py CASE[,CASE] (two processes, two cases each,
take about 20 s). Then v5_median_grid.py."""
import os, sys, json, time
O = os.path.dirname(os.path.abspath(__file__))
os.environ["V5_SCRATCH"] = O + "/v5scratch"; os.makedirs(O + "/v5scratch", exist_ok=True)
os.environ["NUMBA_CACHE_DIR"] = O + "/v5scratch/numba_cache"
R2 = "/home/claude/cfb-wainwright-fire-archive-2"
os.environ["CFB_ARCHIVE_1"] = "/home/claude/cfb-wainwright-fire-archive"
sys.path.insert(0, R2 + "/scripts/rubicon_ridge/model_v5")
import numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
import v5_common as C
REPO_OUT = C.OUT
C.OUT = O + "/v5scratch"; C.SCRATCH = O + "/v5scratch"
from v5_inputs import sample
import geopandas as gpd
from pyproj import Transformer

cases = sys.argv[1].split(",")
S = C.load_v4_setup(); H, W, RES = S["H"], S["W"], S["RES"]; X, Y = C.grid_xy(S)
geom, scar_m = C.load_scar(C.DEFAULT_SCAR, S, "largest")
assert int(scar_m.sum()) == 6750
sb = gpd.read_file(C.SPRING_BURN).to_crs(S["CRS"]); sb_m = S["rast"](sb.union_all())
FP = json.load(open(REPO_OUT + "/front_plane.json")); gx, gy = FP["fit"]["gradient_min_per_km"]
fx, fy = Transformer.from_crs(4326, S["CRS"], always_xy=True).transform(*C.FIRE_LONLAT)
dfront = gx * (X - fx) / 1000 + gy * (Y - fy) / 1000
dy, dx, dist, brg = C.offsets(RES)
runs = pd.read_csv(REPO_OUT + "/ensemble_runs.csv"); runs = runs[runs.weighted & runs.case.isin(cases)]
FR = {"15:27": 27.0, "15:29": 29.0, "15:31": 31.0}
names, grids, t0 = [], [], time.time()
for (case, flab), blk in runs.groupby(["case", "front"]):
    with rasterio.open(REPO_OUT + f"/prefront/arrival_prefront_{case}.tif") as r:
        a = np.full((H, W), np.nan, "float32")
        reproject(r.read(1), a, src_transform=r.transform, src_crs=r.crs, src_nodata=np.nan, dst_transform=S["T"], dst_crs=S["CRS"], dst_nodata=np.nan, resampling=Resampling.nearest)
    pre_min = np.where(np.isfinite(a), (a.astype("float64") - 15) * 60, np.inf)
    tfire = FR[flab]; tf = tfire + dfront
    seed = np.where(pre_min <= tf, pre_min, np.inf); seed = np.where(S["hull_m"], np.minimum(seed, -11.0), seed)
    start = np.isfinite(seed)
    fuel_scar = (S["rast"](geom.buffer(150)) | start) & ~S["water"]
    for _, row in blk.iterrows():
        z = sample(row.schedule, tfire, row.fuel, float(row.ros_mult), row["shape"])
        rfac = np.where(sb_m, float(row.spring_factor), 1.0)
        As = C.spread_v5(seed, tf, fuel_scar, rfac, z["head"], z["ros"], z["ecc"], tfire, 210.0, dy, dx, dist, brg)
        names.append(row.run); grids.append(np.clip(As, -999, 400).astype("float32"))
    print(case, flab, len(names), round(time.time() - t0), "s", flush=True)
np.savez_compressed(O + f"/grids_{'+'.join(cases)}.npz", runs=np.array(names), grids=np.stack(grids).astype("float16") if False else np.stack(grids))
meta = dict(minx=S["minx"], maxy=S["maxy"], res=RES, H=H, W=W, crs=str(S["CRS"]))
json.dump(meta, open(O + "/grid_meta.json", "w"))
print("done", round(time.time() - t0), "s")
