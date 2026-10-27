"""Rubicon Ridge v5, step 2: the four v4 pre-front cases, run on to 15:35 MDT (v4 stopped at 15:16, the front's time at CZWN).
Usage: python3 v5_progression.py ARCHIVE_1 [ARCHIVE_2]
(ARCHIVE_1 = checkout of cfb-wainwright-fire-archive with the v4 code and inputs; outputs go to ARCHIVE_2, default: this repository)
Adapted from scripts/authentic_stoic_alternative_hypothesis/model_v4/progression.py (same part 8b model, same ISI fix, same cases
from model_v4/cases.json and start line from cases_line.json). One change: after 15:00 MDT the pre-front weather is held at the
15:00 observation (CZWN SW wind, about 30 km/h from 234-236 deg). v4 interpolated linearly between the 15:00 and 16:00 hourly
records, so its 15:00-15:16 weather was already drifting toward the post-front 16:00 values; with the front now measured at the
fire at about 15:29, that interpolation would leak post-front wind into the pre-front fire.
Writes data/rubicon-ridge/model_v5/prefront/arrival_prefront_<case>.tif (hours since midnight Oct. 17) and prefront_progression.csv."""
import sys, os, json
import numpy as np, pandas as pd, rasterio
REPO = sys.argv[1]
OUTREPO = sys.argv[2] if len(sys.argv) > 2 else os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../.."))
src = open(os.path.join(REPO, "scripts/authentic_stoic_alternative_hypothesis/part8/ignition_backcalc.py")).read().split("# candidates:")[0]
src = src.replace("return 0.208 * fw * sf", "return fw * sf"); assert "return fw * sf" in src            # v4 ISI fix
src = src.replace("lb = 1.1 * max(ws, 1) ** 0.464", "lb = max(1.0, LBM * 1.1 * max(ws, 1) ** 0.464)")
src = src.replace("def wx(t):\n", "def wx(t):\n    t = min(t, 15.0)   # v5: hold the 15:00 pre-front weather until the front arrives\n"); assert "t = min(t, 15.0)" in src
LBM = 1.0
sys.argv = ["x", REPO, "40", "240"]; exec(src)
V4 = os.path.join(REPO, "data/authentic-stoic-stress-test/model_v4")
OUTD = os.path.join(OUTREPO, "data/rubicon-ridge/model_v5/prefront"); os.makedirs(OUTD, exist_ok=True)
line = set()
LINE = json.load(open(os.path.join(V4, "cases_line.json")))
for dd in np.arange(LINE[1], LINE[2] + 1e-4, RES / 2000):
    lo, la, _ = G.fwd(R11[0], R11[1], LINE[0], dd * 1000); i, j = cell(*tf(lo, la))
    if fuel[i, j]: line.add((i, j))
CASES = [tuple(c) for c in json.load(open(os.path.join(V4, "cases.json")))]
TEND = 15 + 35 / 60
rows = []
for name, mult, lbm, t0 in CASES:
    MORNING_MULT, LBM = mult, lbm
    A = run(0, 0, t0, starts=sorted(line), tend=TEND)
    with rasterio.open(os.path.join(OUTD, f"arrival_prefront_{name}.tif"), "w", driver="GTiff", height=H, width=W, count=1, dtype="float32",
                       crs=CRS, transform=T, nodata=np.nan, compress="deflate") as o:
        o.write(np.where(np.isfinite(A) & (A <= TEND), A, np.nan).astype("float32"), 1)
    r = dict(case=name, ros_multiplier=mult, lb_multiplier=lbm, ignition=f"{int(t0):02d}:{int(round((t0 % 1) * 60)):02d}", K_calibration=round(K, 4))
    for lab, hh in (("1217", TL), ("1449", TV), ("1500", 15.0), ("1510", 15 + 10 / 60), ("1516", 15 + 16 / 60), ("1527", 15 + 27 / 60), ("1529", 15 + 29 / 60), ("1531", 15 + 31 / 60)):
        r[f"area_km2_{lab}"] = round(float((A <= hh).sum()) * RES * RES / 1e6, 3)
    # pre-front growth rate around the GOES 15:10 scan (km2/h), 15:05-15:15
    r["growth_km2_per_h_1505_1515"] = round(float(((A > 15 + 5 / 60) & (A <= 15 + 15 / 60)).sum()) * RES * RES / 1e6 * 6, 3)
    rows.append(r); print(r, flush=True)
pd.DataFrame(rows).to_csv(os.path.join(OUTREPO, "data/rubicon-ridge/model_v5/prefront_progression.csv"), index=False)
H_.round(3).to_csv(os.path.join(OUTD, "prefront_weather_hourly.csv"), index=False)
