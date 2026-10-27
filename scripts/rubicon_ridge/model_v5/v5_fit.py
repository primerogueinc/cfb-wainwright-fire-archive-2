"""Rubicon Ridge v5, step 4: fit weights, arrival-time summaries, GOES checks, figures and the v4 comparison.
Usage:
  python3 v5_fit.py                                   fit to the default scar (scar_patches_revised.geojson, largest patch)
  python3 v5_fit.py --scar PATH --tag NAME [--scar-mode union]
                                                      REFIT to another burn outline (e.g. the Satellogic 70 cm outline) from the
                                                      stored free-burn grids, without re-running the ensemble. Outputs go to
                                                      data/rubicon-ridge/model_v5/scar_<NAME>/.
Reads the per-block checkpoints written by v5_ensemble.py (scratch grids folder: rows_*.csv, free_*.npz, scar_*.npz).

Fit (as v4 fit_post.py): each 'free' run (burns anything but open water) is stopped when its burned area equals the scar's; its
Jaccard overlap J with the scar weights it, w = exp((J - J_max) / 0.05), over the runs with a reduced spring-burn factor
(0.5, 0.75). Runs with factor 1.0 (spring burn treated like any other grass) are a control and get no weight.
Arrival times: from the 'scar-limited' runs (scar + 150 m), the cell containing each reference point and the earliest arrival
within 60 m of it (points 10-20 m outside the mapped edge sit on the edge at 40 m resolution). Free runs are reported as a check.
GOES-18 check (indicative only): modelled fire radiative power from the cells igniting within +/-1 min of each scan's assumed time.
A sensitivity weighting multiplies in a GOES term exp(-0.5 (ln(model/obs) / ln 3)^2) over the 15:10, 15:30 and 15:40 scans.
"""
import os, sys, json, glob, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import v5_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--scar", default=None); ap.add_argument("--scar-mode", default="largest"); ap.add_argument("--tag", default="default")
ap.add_argument("--grids", default=os.path.join(C.SCRATCH, "grids"), help="checkpoint folder written by v5_ensemble.py")
ap.add_argument("--no-figures", action="store_true")
args = ap.parse_args()
OUTD = C.OUT if args.tag == "default" else os.path.join(C.OUT, f"scar_{args.tag}"); os.makedirs(OUTD, exist_ok=True)
SCALE = 0.05

# ---------------------------------------------------------------- load checkpoints
files = sorted(glob.glob(os.path.join(args.grids, "rows_*.csv")))
df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
cfg = json.load(open(os.path.join(C.OUT, "config.json")))
H, W, RES = cfg["grid"]["H"], cfg["grid"]["W"], cfg["grid"]["res_m"]
FRONT_MIN = {"15:27": 27.0, "15:29": 29.0, "15:31": 31.0}
df["front_min"] = df.front.map(FRONT_MIN)
factors = ["case", "front", "schedule", "fuel", "ros_mult", "spring_factor", "shape"]
design = {f: sorted(df[f].unique().tolist(), key=str) for f in factors}
full = len(cfg["cases"]) * len(cfg["front_times_at_fire"]) * len(cfg["schedules"]) * len(cfg["fuels"]) * len(cfg["ros_multipliers"]) * 3 * 2
print("runs", len(df), "blocks", len(files), "design", {k: len(v) for k, v in design.items()}, "complete" if full == len(df) else "PARTIAL")

S = None
def setup():
    global S
    if S is None: S = C.load_v4_setup()
    return S

# ---------------------------------------------------------------- refit to another outline (from stored free grids)
scar_m = np.load(os.path.join(args.grids, "scar_m.npy")) if os.path.exists(os.path.join(args.grids, "scar_m.npy")) else None
if args.scar:
    S = setup(); geom, scar_m = C.load_scar(args.scar, S, args.scar_mode); A = int(scar_m.sum())
    sb_m = np.load(os.path.join(args.grids, "spring_m.npy"))
    J, stop, sbs = {}, {}, {}
    for f in sorted(glob.glob(os.path.join(args.grids, "free_*.npz"))):
        z = np.load(f)
        for k in z.files:
            a = z[k].astype("float32"); fin = np.sort(a[a < 399].ravel())
            ts = float(fin[A - 1]) if len(fin) >= A else np.inf
            b = a <= ts if np.isfinite(ts) else a < 399
            J[k.replace("__", "|")] = float((b & scar_m).sum() / (b | scar_m).sum()); stop[k.replace("__", "|")] = ts
            sbs[k.replace("__", "|")] = float((b & sb_m).sum() / sb_m.sum())
    df["jaccard_default_scar"] = df.jaccard
    df["jaccard"] = df.run.map(J).round(4); df["stop_min"] = df.run.map(stop).round(1); df["stop_time"] = df.stop_min.map(C.hm)
    df["spring_burn_burned_share"] = df.run.map(sbs).round(3)
    for c in [c for c in df.columns if c.startswith("arr_free|")]:
        df[c] = np.nan      # free-run arrivals depend on the stop time; recomputed below from the grids
    scar_info = dict(path=os.path.relpath(args.scar, C.REPO), mode=args.scar_mode, cells=A, km2=round(A * RES * RES / 1e6, 2))
else:
    scar_info = cfg["scar"]

# ---------------------------------------------------------------- weights
mx = df.loc[df.weighted, "jaccard"].max()
df["weight"] = np.where(df.weighted, np.exp((df.jaccard - mx) / SCALE), 0.0); df["weight"] /= df.weight.sum()
GOES_OBS = cfg["frp_check"]["goes_obs_mw"]
lr = sum(((np.log(np.maximum(df[f"frp_model_mw_{k}"], 1.0) / v)) / np.log(3)) ** 2 for k, v in GOES_OBS.items())
df["goes_term"] = np.exp(-0.5 * lr)
df["weight_goes"] = df.weight * df.goes_term; df["weight_goes"] /= df.weight_goes.sum()
mx_all = df.jaccard.max()      # sensitivity: the spring-burn factor 1.0 control runs weighted too
df["weight_incl_control"] = np.exp((df.jaccard - mx_all) / SCALE); df["weight_incl_control"] /= df.weight_incl_control.sum()

# ---------------------------------------------------------------- point arrivals from the stored grids (cell and within 60 m)
S = setup()
PTS = cfg["reference_points"]; PIDX = {k: C.cell_of(S, *v)[:2] for k, v in PTS.items()}
X, Y = C.grid_xy(S)
NEAR = {}
for k, (i, j) in PIDX.items():
    x0, y0 = X[i, j], Y[i, j]; m = np.hypot(X - x0, Y - y0) <= 60.0 + 1e-6
    ii, jj = np.where(m); NEAR[k] = (ii, jj)
near_rows = {}
for f in sorted(glob.glob(os.path.join(args.grids, "scar_*.npz"))):
    z = np.load(f)
    for k in z.files:
        a = z[k].astype("float32"); a[a >= 399] = np.nan
        near_rows[k.replace("__", "|")] = {p: (float(np.nanmin(a[ii, jj])) if np.isfinite(a[ii, jj]).any() else np.nan) for p, (ii, jj) in NEAR.items()}
for p in PTS:
    df[f"arr_scar60|{p}"] = df.run.map(lambda r: near_rows.get(r, {}).get(p, np.nan)).round(2)
if args.scar:   # free-run arrivals at the new stop times
    for f in sorted(glob.glob(os.path.join(args.grids, "free_*.npz"))):
        z = np.load(f)
        for k in z.files:
            run = k.replace("__", "|"); a = z[k].astype("float32"); ts = stop[run]
            for p, (i, j) in PIDX.items():
                v = float(a[i, j]); df.loc[df.run == run, f"arr_free|{p}"] = round(v, 2) if v < 399 and v <= ts else np.nan

df.to_csv(os.path.join(OUTD, "ensemble_runs.csv"), index=False)
wcols = ["run"] + factors + ["weighted", "jaccard", "stop_time", "head_ros_kmh_mean_0_30", "goes_term", "weight", "weight_goes", "weight_incl_control"]
df.sort_values("weight", ascending=False)[wcols].to_csv(os.path.join(OUTD, "fit_weights.csv"), index=False)

# ---------------------------------------------------------------- summaries
def q3(x, w):
    return dict(median=C.wquant(x, w, 0.5), p10=C.wquant(x, w, 0.1), p90=C.wquant(x, w, 0.9))
def rnd(d, n=1): return {k: (round(v, n) if isinstance(v, float) else v) for k, v in d.items()}
def clock(d): return {k: C.hm(v) for k, v in d.items()}
WIN = [("before 3:30", -99, 30), ("3:30-3:45", 30, 45), ("3:45-4:00", 45, 60), ("after 4:00", 60, 999)]
def point_summary(col, sub, wcol):
    x = sub[col].values; w = sub[wcol].values if wcol else np.ones(len(sub))
    ok = np.isfinite(x)
    if not ok.any(): return dict(reached_share=0.0)
    rel = x - sub.front_min.values
    out = dict(reached_share=round(float(w[ok].sum() / w.sum()), 3),
               clock=clock(q3(x, w)), minutes_after_front_at_fire=rnd(q3(rel, w)),
               minutes_after_1500=rnd(q3(x, w)),
               share_by_window={lab: round(float(w[ok & (x > a) & (x <= b)].sum() / w[ok].sum()), 3) for lab, a, b in WIN})
    return out
wsub = df[df.weighted]; csub = df[~df.weighted]
POINTS_MAIN = ["Range 2 (RSO)", "Range 15 (RSO)", "Road head"]
summary = dict(model="Rubicon Ridge, working model v5", scar=scar_info, n_runs=int(len(df)), n_weighted=int(df.weighted.sum()),
               design={k: [str(x) for x in v] for k, v in design.items()}, complete_design=bool(full == len(df)),
               weighting=f"w = exp((J - J_max)/{SCALE}) over runs with spring-burn factor 0.5 or 0.75", J_max=round(float(mx), 4),
               note="Clock times MDT, Oct. 17, 2026. 'cell' = the 40 m cell containing the point; 'within 60 m' = earliest arrival within 60 m; "
                    "'free' = unconstrained runs stopped at the scar's area (point counted only if it burns). Arrival distributions are over the "
                    "weighted runs; 'equal' = all runs with spring-burn factor 0.5/0.75 equally weighted.")
summary["points"] = {}
for p in PTS:
    summary["points"][p] = dict(
        scar_limited_cell_weighted=point_summary(f"arr_scar|{p}", wsub, "weight"),
        scar_limited_within60m_weighted=point_summary(f"arr_scar60|{p}", wsub, "weight"),
        scar_limited_cell_equal=point_summary(f"arr_scar|{p}", wsub, None),
        scar_limited_cell_goes_weighted=point_summary(f"arr_scar|{p}", wsub, "weight_goes"),
        free_cell_weighted=point_summary(f"arr_free|{p}", wsub, "weight"),
        by_front_weighted={fr: point_summary(f"arr_scar|{p}", wsub[wsub.front == fr], "weight") for fr in FRONT_MIN},
        control_spring_factor_1_equal=point_summary(f"arr_scar|{p}", csub, None),
        all_runs_incl_control_weighted=point_summary(f"arr_scar|{p}", df, "weight_incl_control"),
        by_ros_mult_weighted={str(m): point_summary(f"arr_scar|{p}", wsub[wsub.ros_mult == m], "weight") for m in sorted(wsub.ros_mult.unique())},
        head_ros_ge_11kmh_weighted=point_summary(f"arr_scar|{p}", wsub[wsub.head_ros_kmh_mean_0_30 >= 11], "weight"))
# marginal weights
summary["weight_by_factor"] = {f: {str(k): round(float(v), 3) for k, v in wsub.groupby(f).weight.sum().items()} for f in factors if f != "spring_factor"}
summary["weight_by_factor"]["spring_factor"] = {str(k): round(float(v), 3) for k, v in wsub.groupby("spring_factor").weight.sum().items()}
summary["weight_by_factor_incl_control"] = {f: {str(k): round(float(v), 3) for k, v in df.groupby(f).weight_incl_control.sum().items()} for f in factors}
summary["weight_by_factor_goes"] = {f: {str(k): round(float(v), 3) for k, v in wsub.groupby(f).weight_goes.sum().items()} for f in factors}
summary["jaccard_by_factor_max"] = {f: {str(k): round(float(v), 3) for k, v in df.groupby(f).jaccard.max().items()} for f in factors}
summary["implied_post_front_ros_kmh"] = dict(
    mean_first_30_min_weighted=rnd(q3(wsub.head_ros_kmh_mean_0_30, wsub.weight), 2),
    at_front_weighted=rnd(q3(wsub.head_ros_kmh_at_front, wsub.weight), 2),
    mean_first_60_min_weighted=rnd(q3(wsub.head_ros_kmh_mean_0_60, wsub.weight), 2),
    ensemble_range_first_30_min=[round(float(df.head_ros_kmh_mean_0_30.min()), 2), round(float(df.head_ros_kmh_mean_0_30.max()), 2)],
    reference=dict(fbp_O1a_max_kmh=11.4, fbp_O1b_approx_kmh=15, independent_expert_expected_kmh=12.7),
    share_weight_ros_ge_11=round(float(wsub.weight[wsub.head_ros_kmh_mean_0_30 >= 11].sum()), 3))
best = df[df.weighted].sort_values("jaccard", ascending=False).head(10)
summary["best_runs"] = best[["run", "jaccard", "weight", "stop_time", "head_ros_kmh_mean_0_30", "spring_burn_burned_share"] +
                            [f"arr_scar|{p}" for p in POINTS_MAIN]].round(4).to_dict("records")
summary["observed_spring_burn_share_burned"] = cfg["spring_burn"]["share_burned_in_fit_scar"] if not args.scar else round(float((np.load(os.path.join(args.grids, "spring_m.npy")) & scar_m).sum() / np.load(os.path.join(args.grids, "spring_m.npy")).sum()), 3)
summary["spring_burn_share_burned_weighted"] = rnd(q3(wsub.spring_burn_burned_share, wsub.weight), 3)
summary["control_spring_factor_1"] = dict(J_max=round(float(csub.jaccard.max()), 4), J_median=round(float(csub.jaccard.median()), 4),
                                          J_max_weighted_runs=round(float(mx), 4), J_median_weighted_runs=round(float(wsub.jaccard.median()), 4))
# GOES check
g = {}
for k, v in GOES_OBS.items():
    col = f"frp_model_mw_{k}"
    g[k] = dict(observed_mw=v, model_weighted=rnd(q3(wsub[col], wsub.weight), 0),
                model_by_front_median={fr: (round(float(wsub.loc[wsub.front == fr, col].median())) if (wsub.front == fr).any() else None) for fr in FRONT_MIN})
pre = df.groupby(["case", "front"]).agg(prefront_growth_km2_1510_to_front=("prefront_growth_km2_1510_to_front", "first"),
                                         frp_model_mw_1510=("frp_model_mw_15:10", "median"), frp_1500_1520_mean=("frp_model_mw_1500_1520_mean", "median")).reset_index()
g["prefront_by_case_front"] = pre.round(3).to_dict("records")
g["no_large_growth_before_front"] = bool((pre.frp_model_mw_1510 <= 3 * GOES_OBS["15:10"]).all())
g["note"] = ("Indicative: FRP from cells igniting within +/-1 min of the scan's assumed time (15:11, 15:31, 15:41), 0.35 kg/m2, 18 MJ/kg, "
             "radiative fraction 0.15. GOES-18 reads 0 MW at 15:20 (cloud/plume obscuration), so single scans carry large uncertainty.")
summary["goes_check"] = g
json.dump(summary, open(os.path.join(OUTD, "arrival_summary.json"), "w"), indent=1, default=lambda o: None if o is None else float(o))

# ---------------------------------------------------------------- v4 vs v5 table
V4 = os.path.join(C.V4REPO, "data/authentic-stoic-stress-test/model_v4")
e4 = json.load(open(os.path.join(V4, "ensemble_summary.json"))); o4 = json.load(open(os.path.join(V4, "overtaken_time.json"))); o4b = json.load(open(os.path.join(V4, "overtaken_time_f325.json")))
P = summary["points"]
def fmt(ps):
    if ps.get("reached_share", 0) == 0: return "not reached"
    c = ps["clock"]; m = ps["minutes_after_front_at_fire"]
    return f"{c['median']} ({c['p10']}-{c['p90']}); +{m['median']} min ({m['p10']}-{m['p90']})"
v4c = lambda x: "15" + x[1:] if x.startswith("3:") else x       # v4 files use 12-hour clock (3:21 = 15:21)
rows = [
    ("Front at the fire", "15:16 (CZWN station) or 15:25 (allowance)", "15:27 / 15:29 / 15:31 (moving line fitted to 9 stations/sensors; 15:28 +/- 0.5 best)"),
    ("Front geometry", "same time everywhere", f"moving line toward {cfg['front_plane']['toward_deg']} deg at {cfg['front_plane']['speed_kmh']} km/h"),
    ("Post-front wind", "5 schedules incl. 'authentic stoic' slow veer; fixed 306/315 deg", "5 schedules from SUSTAINED 10 m wind: CZWN observed veer 290-307 deg, fixed 300/310/315 deg, RAP/HRDPS mean"),
    ("Post-front spread", "head ROS fixed 4.4 / 8 / 12 km/h", "FBP O-1a or O-1b from ISI (hourly FFMC) x 0.5/1.0/1.25 (about 3.4-12.9 km/h)"),
    ("April 20 spring burn", "same as other grass", "spread factor 0.5 or 0.75 (1.0 as unweighted control)"),
    ("Runs", str(e4["n_fires"]), str(len(df))),
    ("Best overlap with scar (Jaccard)", "0.51-0.55", f"{mx:.3f}"),
    ("Range 2, v4 OSM point (weighted median, 10-90%)", f"{v4c(o4['points_weighted_share_by_window']['Range 2 point']['median_arrival_clock'])} (front 15:16) / {v4c(o4b['points_weighted_share_by_window']['Range 2 point']['median_arrival_clock'])} (front 15:25); +{e4['overall_scar_weighted']['arrival_range2_min']['median']} min ({e4['overall_scar_weighted']['arrival_range2_min']['p10']}-{e4['overall_scar_weighted']['arrival_range2_min']['p90']})",
     fmt(P["Range 2 point (OSM, v4)"]["scar_limited_cell_weighted"])),
    ("Range 2, official RSO grid (weighted)", "not modelled", fmt(P["Range 2 (RSO)"]["scar_limited_cell_weighted"])),
    ("Range 15 (weighted)", "not modelled", fmt(P["Range 15 (RSO)"]["scar_limited_cell_weighted"])),
    ("Road head (weighted)", f"{v4c(o4['points_weighted_share_by_window']['Road head']['median_arrival_clock'])} (front 15:16) / {v4c(o4b['points_weighted_share_by_window']['Road head']['median_arrival_clock'])} (front 15:25); +{e4['overall_scar_weighted']['arrival_roadhead_min']['median']} min ({e4['overall_scar_weighted']['arrival_roadhead_min']['p10']}-{e4['overall_scar_weighted']['arrival_roadhead_min']['p90']})",
     fmt(P["Road head"]["scar_limited_cell_weighted"])),
    ("Share of weight after 15:45 at road head", f"{o4['points_weighted_share_by_window']['Road head']['3:45-4:00'] + o4['points_weighted_share_by_window']['Road head']['after 4:00']:.2f} (front 15:16) / {o4b['points_weighted_share_by_window']['Road head']['3:45-4:00'] + o4b['points_weighted_share_by_window']['Road head']['after 4:00']:.2f} (front 15:25)",
     f"{P['Road head']['scar_limited_cell_weighted']['share_by_window']['3:45-4:00'] + P['Road head']['scar_limited_cell_weighted']['share_by_window']['after 4:00']:.2f}"),
]
cmp = pd.DataFrame(rows, columns=["item", "v4 (Oct. 25)", "v5 Rubicon Ridge"]); cmp.to_csv(os.path.join(OUTD, "v4_vs_v5.csv"), index=False)
print(cmp.to_string(index=False))
print(json.dumps(summary["weight_by_factor"], indent=0)); print(json.dumps(summary["implied_post_front_ros_kmh"]))
for p in POINTS_MAIN + ["Range 2 point (OSM, v4)"]:
    print(p, json.dumps(P[p]["scar_limited_cell_weighted"]), "| 60m", json.dumps(P[p]["scar_limited_within60m_weighted"]["clock"]), "| free reached", P[p]["free_cell_weighted"].get("reached_share"))
print(json.dumps(g, default=float)[:1500])
if args.no_figures: sys.exit(0)

# ---------------------------------------------------------------- figures
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
INK, MUTED, GRID = "#1f2328", "#59636e", "#d0d7de"
SER = ["#2f6db5", "#c2552d", "#3d8b5a"]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "font.family": "DejaVu Sans"})
sb_m = np.load(os.path.join(args.grids, "spring_m.npy")); water = np.load(os.path.join(args.grids, "water.npy"))
ext = [S["minx"], S["minx"] + W * RES, S["maxy"] - H * RES, S["maxy"]]
def grid_of(run, kind):
    key = "_".join([run.split("|")[0], run.split("|")[1].replace(":", "")])
    return np.load(os.path.join(args.grids, f"{kind}_{key}.npz"))[run.replace("|", "__")].astype("float32")
# best run of each of the three best distinct configurations (wind schedule, shape, spring factor): near-duplicates differ only in
# spread rate or front time, which the area-matched fit cannot tell apart
top = df[df.weighted].sort_values("jaccard", ascending=False).drop_duplicates(["schedule", "shape", "spring_factor"]).head(3)
fig, axs = plt.subplots(1, 3, figsize=(13, 4.9), constrained_layout=True)
pts_xy = {p: C.cell_of(S, *PTS[p])[2:] for p in POINTS_MAIN}
for ax, (_, r) in zip(axs, top.iterrows()):
    a = grid_of(r.run, "free"); b = a <= r.stop_min
    ax.imshow(np.where(water, 1, np.nan), extent=ext, cmap="Blues", vmin=0, vmax=3, interpolation="nearest")
    ax.imshow(np.where(b, 1, np.nan), extent=ext, cmap="Oranges", vmin=0, vmax=2.2, alpha=0.75, interpolation="nearest")
    ax.contour(X, Y, scar_m.astype(float), levels=[0.5], colors=INK, linewidths=1.0)
    ax.contour(X, Y, sb_m.astype(float), levels=[0.5], colors="#3d8b5a", linewidths=0.8, linestyles="--")
    ax.contour(X, Y, np.where(b, a, np.nan), levels=[r.front_min + 10, r.front_min + 30, r.front_min + 60], colors=MUTED, linewidths=0.5)
    for p, (x, y) in pts_xy.items():
        ax.plot(x, y, "o", ms=4, mfc="white", mec=INK, mew=0.9)
    ax.set_xlim(ext[0] + 1500, ext[1] - 1500); ax.set_ylim(ext[2] + 1500, ext[3] - 1500); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"J = {r.jaccard:.3f}  ({r.front} front, {r.schedule.split('_')[1]}, {r.fuel} x{r.ros_mult}, spring x{r.spring_factor}, {r["shape"]})\n"
                 f"stopped at {r.stop_time}; head ROS {r.head_ros_kmh_mean_0_30:.1f} km/h (first 30 min)", fontsize=8, color=INK, loc="left")
axs[0].legend(handles=[Line2D([], [], color=INK, lw=1, label=f"burn scar ({os.path.basename(scar_info['path'])})"),
                       Line2D([], [], color="#3d8b5a", lw=0.8, ls="--", label="April 20, 2026 planned burn"),
                       Line2D([], [], color="#e8873a", lw=6, alpha=0.75, label="modelled burn (stopped at scar area)"),
                       Line2D([], [], color=MUTED, lw=0.5, label="isochrones: front +10 / +30 / +60 min"),
                       Line2D([], [], ls="", marker="o", mfc="white", mec=INK, label="Range 2, Range 15, road head")],
              loc="lower left", fontsize=7, frameon=False)
fig.suptitle("Rubicon Ridge v5: best-fitting post-front runs (three distinct wind / shape / spring-burn set-ups) against the burn scar", x=0.01, ha="left", color=INK, fontsize=11)
fig.savefig(os.path.join(OUTD, "fig_best_runs_vs_scar.png"), dpi=150); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True, constrained_layout=True)
bins = np.arange(20, 66, 1.0)
for ax, p in zip(axs, POINTS_MAIN):
    for c, fr in zip(SER, FRONT_MIN):
        s = wsub[wsub.front == fr]; x = s[f"arr_scar|{p}"].values; ok = np.isfinite(x)
        ax.hist(x[ok], bins=bins, weights=s.weight.values[ok], histtype="step", lw=1.6, color=c, label=f"front at fire {fr}")
    ps = P[p]["scar_limited_cell_weighted"]
    for qk, ls in (("p10", ":"), ("median", "-"), ("p90", ":")):
        v = ps["minutes_after_1500"][qk]; ax.axvline(v, color=INK, lw=0.9, ls=ls)
    ax.set_title(f"{p}: median {ps['clock']['median']} (10-90%: {ps['clock']['p10']}-{ps['clock']['p90']})", fontsize=9, loc="left", color=INK)
    ax.set_xticks(range(20, 66, 5)); ax.set_xticklabels([f"3:{m:02d}" if m < 60 else f"4:{m - 60:02d}" for m in range(20, 66, 5)])
    ax.grid(axis="y", color=GRID, lw=0.5); ax.spines[["top", "right"]].set_visible(False); ax.set_xlabel("fire arrival, p.m. MDT")
axs[0].set_ylabel("share of fit weight per minute"); axs[0].legend(fontsize=7, frameon=False)
fig.suptitle("Rubicon Ridge v5: when the modelled fire reaches each point (weighted by fit to the scar; scar-limited runs)", x=0.01, ha="left", color=INK, fontsize=11)
fig.savefig(os.path.join(OUTD, "fig_arrival_distributions.png"), dpi=150); plt.close(fig)

import textwrap
wrapw = [26, 46, 62]
cells = [["\n".join(textwrap.wrap(str(v), w)) for v, w in zip(row, wrapw)] for row in cmp.values]
nl = [max(c.count("\n") + 1 for c in row) for row in cells]
fig, ax = plt.subplots(figsize=(13, 0.3 * sum(nl) + 1.2)); ax.axis("off")
t = ax.table(cellText=cells, colLabels=cmp.columns, loc="upper center", cellLoc="left", colWidths=[0.2, 0.34, 0.46])
t.auto_set_font_size(False); t.set_fontsize(8)
for (i, j), c in t.get_celld().items():
    c.set_edgecolor(GRID); c.set_height((nl[i - 1] if i > 0 else 1) / (sum(nl) + 1.5))
    if i == 0: c.set_text_props(weight="bold", color=INK); c.set_facecolor("#f3f5f7")
ax.set_title("Working model v4 vs Rubicon Ridge v5", loc="left", color=INK, fontsize=11)
fig.savefig(os.path.join(OUTD, "fig_v4_vs_v5.png"), dpi=150, bbox_inches="tight"); plt.close(fig)
print("figures written to", OUTD)
