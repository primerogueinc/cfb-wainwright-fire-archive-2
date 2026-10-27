"""Rubicon Ridge v5, step 3: the post-front ensemble, its fit to the burn scar and the GOES checks.
Usage: python3 v5_ensemble.py [--scar PATH] [--scar-mode largest|union] [--tag NAME] [--edge DIR]
  --scar       burn scar used for fitting (default data/infrastructure/scar_patches_revised.geojson, as v4). Swap in a sharper
               outline (e.g. the Satellogic 70 cm outline; faster: v5_fit.py --scar, which refits from the stored grids) here; the grid, starting fire and everything else stay fixed.
  --scar-mode  'largest' polygon (v4 behaviour, default) or 'union' of all polygons.
  --tag        output sub-folder name (default 'default'); outputs go to data/rubicon-ridge/model_v5/ (tag 'default') or
               data/rubicon-ridge/model_v5/scar_<tag>/.
  --edge       folder with v4 edge_check.py grids (grid_meta.npy, S2_20_m_dNBR_2026_10_24.npy) for the road-head pocket test; optional.

Each run (one combination of the factors in config.json):
  starting state: the v5 pre-front case (v5_progression.py) burns until the front reaches each cell; the outline of the 12:17
                  Landsat and 2:49 p.m. VIIRS detections is always burned (as v4).
  front:          moving line, t_front(x, y) = t_fire + gradient . (x - fire) (v5_inputs.py, front_plane.json).
  after the front: wind schedule W1-W5, FBP O-1a or O-1b head spread from ISI (hourly FFMC, sustained wind), times a spread
                  multiplier (0.5 / 1.0 / 1.25), shape steady (FBP L:B from the wind) or gusty (L:B 3), and the April 20 spring
                  burn at a reduced spread factor (0.5 or 0.75; 1.0 as an unweighted control).
  checkpoints:    each (case, front time) block of runs is saved to the scratch grids folder (rows_<case>_<front>.csv plus
                  compressed float16 arrival grids); blocks already on disk are skipped, so a restart resumes. Weights, summaries and
                  figures are made afterwards by v5_fit.py.
  two versions:   'free' (anything but open water burns; stopped when its area equals the scar's: Jaccard overlap = fit) and
                  'scar-limited' (scar + 150 m, as v4's timing runs: arrival times at the reference points and the GOES check).
"""
import os, sys, json, math, argparse, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
import v5_common as C
from v5_inputs import sample, SCHED

ap = argparse.ArgumentParser()
ap.add_argument("--scar", default=C.DEFAULT_SCAR); ap.add_argument("--scar-mode", default="largest")
ap.add_argument("--tag", default="default"); ap.add_argument("--edge", default=os.path.join(os.path.dirname(C.SCRATCH), "v4/edge"))
ap.add_argument("--cases", default="", help="comma-separated subset of v4 cases (default all); used to run cases in separate processes")
ap.add_argument("--schedules", default="", help="comma-separated subset of wind schedules (reduced design)")
args = ap.parse_args()
OUTD = C.OUT if args.tag == "default" else os.path.join(C.OUT, f"scar_{args.tag}"); os.makedirs(OUTD, exist_ok=True)
GRIDS = os.path.join(C.SCRATCH, "grids" if args.tag == "default" else f"grids_{args.tag}"); os.makedirs(GRIDS, exist_ok=True)

S = C.load_v4_setup(); H, W, RES = S["H"], S["W"], S["RES"]
scar_geom, scar_m = C.load_scar(args.scar, S, args.scar_mode)
A_SCAR = int(scar_m.sum())
water = S["water"]; X, Y = C.grid_xy(S)

# spring burn (April 20-21, 2026) on the grid
import geopandas as gpd
sb = gpd.read_file(C.SPRING_BURN).to_crs(S["CRS"]); sb_m = S["rast"](sb.union_all())
obs_sb_share = float((sb_m & scar_m).sum() / max(sb_m.sum(), 1))

# moving front
FP = json.load(open(os.path.join(C.OUT, "front_plane.json")))
gx, gy = FP["fit"]["gradient_min_per_km"]
from pyproj import Transformer
fx, fy = Transformer.from_crs(4326, S["CRS"], always_xy=True).transform(*C.FIRE_LONLAT)
dfront = gx * (X - fx) / 1000 + gy * (Y - fy) / 1000          # minutes relative to the time at the fire centre

PTS = C.ref_points(); PIDX = {k: C.cell_of(S, *v)[:2] for k, v in PTS.items()}
dy, dx, dist, brg = C.offsets(RES)

# optional: v4 road-head pocket and edge agreement (Oct. 24 Sentinel-2 dNBR)
pocket = near = obs_b = obs_ok = None
if os.path.exists(os.path.join(args.edge, "grid_meta.npy")):
    from rasterio.transform import from_origin
    from scipy import ndimage
    ex0, ey0, er, eh, ew = np.load(os.path.join(args.edge, "grid_meta.npy"))
    ob10 = np.load(os.path.join(args.edge, "S2_20_m_dNBR_2026_10_24.npy")).astype("float32")
    obs = np.full((H, W), np.nan, "float32")
    reproject(ob10, obs, src_transform=from_origin(ex0, ey0, er, er), src_crs=S["CRS"], dst_transform=S["T"], dst_crs=S["CRS"], resampling=Resampling.average, dst_nodata=np.nan)
    obs_b = obs >= 0.5; obs_ok = np.isfinite(obs)
    rh = C.cell_of(S, *PTS["Road head"]); r2 = C.cell_of(S, *PTS["Range 2 (RSO)"])
    near = (np.hypot(X - rh[2], Y - rh[3]) <= 600) | (np.hypot(X - r2[2], Y - r2[3]) <= 600)
    lab, _ = ndimage.label(~obs_b & obs_ok & near)
    i, j = rh[:2]
    pocket = lab == (lab[i, j] if lab[i, j] else lab.flat[np.argmin(np.where(lab > 0, np.hypot(X - rh[2], Y - rh[3]), np.inf))])

CASES = [c[0] for c in json.load(open(os.path.join(C.V4REPO, "data/authentic-stoic-stress-test/model_v4/cases.json")))]
FRONTS = {"15:27": 27.0, "15:29": 29.0, "15:31": 31.0}
FUELS = ["O-1a", "O-1b"]; MULTS = [0.5, 1.0, 1.25]; SBF = [0.5, 0.75, 1.0]; SHAPES = ["steady", "gusty"]
LOAD = 0.35; HEAT = 18000.0; CHI = 0.15          # kg/m2 consumed, kJ/kg, radiative fraction (FRP check, indicative)
GOES_T = {"15:10": 11.0, "15:30": 31.0, "15:40": 41.0}   # scan start + ~1 min (fire's line in the full-disk scan; assumption)
GOES_OBS = {"15:10": 76.5, "15:30": 1215.4, "15:40": 336.3}
cell_e = RES * RES * LOAD * HEAT * CHI / 1000.0          # MJ radiated per fully burned cell
def frp(A, t, half=1.0, mask=None):
    m = (A > t - half) & (A <= t + half)
    if mask is not None: m &= mask
    e = (np.where(sb_m, sbf_cur, 1.0) * m).sum() * cell_e
    return float(e / (2 * half * 60))                    # MW

config = dict(model="Rubicon Ridge, working model v5", grid=dict(res_m=RES, H=H, W=W, crs=S["CRS"], minx=S["minx"], maxy=S["maxy"]),
              scar=dict(path=os.path.relpath(args.scar, C.REPO), mode=args.scar_mode, cells=A_SCAR, km2=round(A_SCAR * RES * RES / 1e6, 2)),
              cases=CASES, front_times_at_fire=list(FRONTS), front_plane=FP["fit"], schedules={k: v["note"] for k, v in SCHED.items()},
              fuels=FUELS, ros_multipliers=MULTS, spring_burn_factors=dict(weighted=[0.5, 0.75], control_unweighted=[1.0]), shapes=dict(steady="FBP grass L:B = 1.1 WS^0.464", gusty="L:B 3.0"),
              curing_pct=95, weighting="w = exp((J - J_max) / 0.05) over the weighted runs (spring-burn factor 0.5 or 0.75)",
              frp_check=dict(load_kg_m2=LOAD, heat_kJ_kg=HEAT, radiative_fraction=CHI, goes_line_time_assumed=GOES_T, goes_obs_mw=GOES_OBS),
              spring_burn=dict(path=os.path.relpath(C.SPRING_BURN, C.REPO), cells=int(sb_m.sum()), share_burned_in_fit_scar=round(obs_sb_share, 3)),
              reference_points={k: v for k, v in PTS.items()})
json.dump(config, open(os.path.join(OUTD, "config.json"), "w"), indent=1, default=float)

if args.cases: CASES = [c for c in CASES if c in args.cases.split(",")]
if args.schedules: SCHED = {k: v for k, v in SCHED.items() if k in args.schedules.split(",")}
rows = []; t0 = time.time()
for case in CASES:
    with rasterio.open(os.path.join(C.OUT, "prefront", f"arrival_prefront_{case}.tif")) as r:
        a = np.full((H, W), np.nan, "float32")
        reproject(r.read(1), a, src_transform=r.transform, src_crs=r.crs, src_nodata=np.nan, dst_transform=S["T"], dst_crs=S["CRS"], dst_nodata=np.nan, resampling=Resampling.nearest)
    pre_min = np.where(np.isfinite(a), (a.astype("float64") - 15) * 60, np.inf)        # minutes after 15:00 (float64: see spread_v5)
    for flab, tfire in FRONTS.items():
        key = f"{case}_{flab.replace(':', '')}"
        ck = os.path.join(GRIDS, f"rows_{key}.csv")
        if os.path.exists(ck) and os.path.exists(os.path.join(GRIDS, f"scar_{key}.npz")):
            print("skip (checkpoint)", key, flush=True); continue
        rows = []
        tf = tfire + dfront
        seed = np.where(pre_min <= tf, pre_min, np.inf)
        seed = np.where(S["hull_m"], np.minimum(seed, -11.0), seed)   # observed by 2:49 p.m.
        start = np.isfinite(seed)
        fuel_free = ~water; fuel_scar = (S["rast"](scar_geom.buffer(150)) | start) & ~water
        pre_growth = float(((pre_min > 10) & (pre_min <= tf) & ~S["hull_m"]).sum() * RES * RES / 1e6)   # km2 burned 15:10 -> front, pre-front model
        pack_free, pack_scar = {}, {}
        for sch, fuel, mult, sbf, shape in itertools.product(SCHED, FUELS, MULTS, SBF, SHAPES):
            sbf_cur = sbf
            z = sample(sch, tfire, fuel, mult, shape)
            rfac = np.where(sb_m, sbf, 1.0)
            run = f"{case}|{flab}|{sch}|{fuel}|x{mult}|sb{sbf}|{shape}"
            Af = C.spread_v5(seed, tf, fuel_free, rfac, z["head"], z["ros"], z["ecc"], tfire, 240.0, dy, dx, dist, brg)
            As = C.spread_v5(seed, tf, fuel_scar, rfac, z["head"], z["ros"], z["ecc"], tfire, 210.0, dy, dx, dist, brg)
            fin = np.sort(Af[np.isfinite(Af)])
            if len(fin) < A_SCAR: tstop = np.inf; burned = np.isfinite(Af)
            else: tstop = float(fin[A_SCAR - 1]); burned = Af <= tstop
            inter = (burned & scar_m).sum(); uni = (burned | scar_m).sum()
            k30 = int(30 / C.STEP)
            r = dict(run=run, case=case, front=flab, schedule=sch, fuel=fuel, ros_mult=mult, spring_factor=sbf, shape=shape,
                     weighted=sbf < 1.0, jaccard=round(float(inter / uni), 4), stop_time=C.hm(tstop), stop_min=round(tstop, 1),
                     spring_burn_burned_share=round(float((burned & sb_m).sum() / sb_m.sum()), 3),
                     head_ros_kmh_at_front=round(float(z["ros"][0] * 0.06), 2), head_ros_kmh_mean_0_30=round(float(z["ros"][:k30 + 1].mean() * 0.06), 2),
                     head_ros_kmh_mean_0_60=round(float(z["ros"][:int(60 / C.STEP) + 1].mean() * 0.06), 2),
                     wind_kmh_0_30=round(float(z["ws"][:k30 + 1].mean()), 1), isi_0_30=round(float(z["isi"][:k30 + 1].mean()), 1),
                     prefront_growth_km2_1510_to_front=round(pre_growth, 3))
            if pocket is not None:
                r["edge_agree"] = round(float((burned[near & obs_ok] == obs_b[near & obs_ok]).mean()), 3)
                r["pocket_burned"] = round(float(burned[pocket].mean()), 3)
            for k, (i, j) in PIDX.items():
                r[f"arr_scar|{k}"] = round(float(As[i, j]), 2) if np.isfinite(As[i, j]) else np.nan
                r[f"arr_free|{k}"] = round(float(Af[i, j]), 2) if np.isfinite(Af[i, j]) and Af[i, j] <= tstop else np.nan
            for lab, tt in GOES_T.items():
                r[f"frp_model_mw_{lab}"] = round(frp(As, tt))
                for pk, pm in S["foot"].items():
                    r[f"frp_model_mw_{lab}_{pk.split()[0]}"] = round(frp(As, tt, mask=pm))
            r["frp_model_mw_1500_1520_mean"] = round(float(np.mean([frp(As, tt, half=0.5) for tt in np.arange(0.5, 20, 1.0)])))
            for tt, lab in ((31, "1531"), (41, "1541"), (60, "1600")):
                r[f"area_km2_{lab}"] = round(float((As <= tt).sum() * RES * RES / 1e6), 3)
            rows.append(r)
            pack_free[run] = np.clip(Af, -60, 400).astype("float16"); pack_scar[run] = np.clip(As, -60, 400).astype("float16")
        np.savez_compressed(os.path.join(GRIDS, f"free_{key}.npz"), **{k.replace("|", "__"): v for k, v in pack_free.items()})
        np.savez_compressed(os.path.join(GRIDS, f"scar_{key}.npz"), **{k.replace("|", "__"): v for k, v in pack_scar.items()})
        pd.DataFrame(rows).to_csv(ck + ".tmp", index=False); os.replace(ck + ".tmp", ck)     # checkpoint written last
        print(case, flab, len(rows), round(time.time() - t0), "s", flush=True)
np.save(os.path.join(GRIDS, "meta.npy"), np.array([S["minx"], S["maxy"], RES, H, W]))
np.save(os.path.join(GRIDS, "scar_m.npy"), scar_m); np.save(os.path.join(GRIDS, "spring_m.npy"), sb_m); np.save(os.path.join(GRIDS, "water.npy"), water)
print("done; run v5_fit.py --tag", args.tag)
