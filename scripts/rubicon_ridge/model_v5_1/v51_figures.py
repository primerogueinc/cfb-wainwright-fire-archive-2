"""Rubicon Cloud (working model v5.1), step 2: v5 isochrones (vectors) and the figures.
Usage:  python3 v51_figures.py [--gif]      (run v51_air.py first)
  * v5_isochrones_weighted_median.geojson: per-cell weighted median arrival of the v5 scar-limited runs (default scar weights),
    contoured at 15:30-16:00 MDT. Needs the stored v5 grids (V5_SCRATCH/grids); if they are absent the committed file is reused.
  * fig_tracks_map.png, fig_time_alt_dist.png, fig_response_timeline.png, and with --gif anim_tracks_1500-0400.gif.
No imagery of any kind is used. The burn outline is the vector outline in data/rubicon-ridge/satellogic_2026-10-20/.
"""
import os, sys, json, glob, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import v51_common as C
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import geopandas as gpd
from pyproj import Transformer

ap = argparse.ArgumentParser(); ap.add_argument("--gif", action="store_true")
ap.add_argument("--grids", default=os.path.join(C.C5.SCRATCH, "grids")); args = ap.parse_args()
OUT = C.OUT
INK, MUTED, GRID, FIRE_RED = "#1f2328", "#59636e", "#d0d7de", "#b3361b"
OPS = {"STAR": ("#c2552d", "STARS helicopter (call sign STARnn)"), "CWA": ("#2f6db5", "Can West, AHS air ambulance contractor (CWAnnn)"),
       "SLG": ("#3d8b5a", "Sask. Government Air Ambulance (SLGn)"), "": ("#8c6bb1", "operator not identified from call sign")}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "font.family": "DejaVu Sans"})
FOOT = ("Times MDT (UTC-6), Oct. 17-18, 2026. Aircraft: adsb.lol open archive, (c) adsb.lol contributors, ODbL 1.0; only aircraft "
        "that broadcast ADS-B and were heard are shown (military aircraft usually are not).\nOperators are inferred from public call-sign "
        "prefixes only; nothing here says who was aboard. Burn outline: Prime Rogue vectors from Satellogic 70 cm imagery via SkyFi, "
        "Oct. 20, 2026 (no imagery shown). Model: Rubicon Ridge v5. Graphic: Prime Rogue Inc.")
T0 = pd.Timestamp("2026-10-17 15:00"); T1 = pd.Timestamp("2026-10-18 04:30")

def opcode(cs):
    for k in ("STAR", "CWA", "SLG"):
        if isinstance(cs, str) and cs.startswith(k):
            return k
    return ""

# ---------------------------------------------------------------- data
ctx = pd.read_csv(os.path.join(OUT, "context_tracks.csv"), parse_dates=["mdt"])
ctx["alt_ft"] = pd.to_numeric(ctx.alt_baro, errors="coerce")
S = pd.read_csv(C.SUPPLIED, parse_dates=["mdt"])
FL = pd.read_csv(os.path.join(OUT, "flights.csv"), parse_dates=["first_heard_mdt", "last_heard_mdt", "departure_mdt", "arrival_mdt"])
V = pd.read_csv(os.path.join(OUT, "wainwright_visits.csv"), parse_dates=["last_or_first_heard_mdt", "est_on_ground_mdt"])
FR = pd.read_csv(os.path.join(OUT, "inbound_response_flights.csv"), parse_dates=["departure_or_first_heard_mdt", "arrival_mdt"])
AS = pd.read_csv(os.path.join(OUT, "aircraft_summary.csv"))
v5 = json.load(open(os.path.join(C.V5, "arrival_summary.json")))
cfg = json.load(open(os.path.join(C.V5, "config.json")))
burn = gpd.read_file(C.BURN).to_crs(4326)
WIN = {lab: v5["points"][p]["scar_limited_cell_weighted"]["clock"] for p, lab in
       [("Range 2 (RSO)", "Range 2"), ("Range 15 (RSO)", "Range 15"), ("Road head", "road head")]}
hm = lambda s: pd.Timestamp("2026-10-17 " + s)
FRONT = T0 + pd.Timedelta(minutes=C.FRONT_AT_FIRE)
cs_of = {r.registration: (r.callsigns_in_archive if isinstance(r.callsigns_in_archive, str) else "") for r in AS.itertuples()}

# ---------------------------------------------------------------- v5 isochrones (weighted median arrival per cell)
ISO = os.path.join(OUT, "v5_isochrones_weighted_median.geojson")
LEVELS = [30, 35, 40, 45, 50, 60]
if os.path.isdir(args.grids) and glob.glob(os.path.join(args.grids, "scar_*.npz")):
    runs = pd.read_csv(os.path.join(C.V5, "ensemble_runs.csv"), usecols=["run", "weight"])
    wmap = dict(zip(runs.run.str.replace("|", "__", regex=False), runs.weight))
    H, W, RES = cfg["grid"]["H"], cfg["grid"]["W"], cfg["grid"]["res_m"]
    edges = np.arange(-61, 242, 1.0)                     # 1-min bins, minutes after 15:00; last bin = later / never
    hist = np.zeros((len(edges), H * W)); tot = 0.0
    for f in sorted(glob.glob(os.path.join(args.grids, "scar_*.npz"))):
        z = np.load(f)
        for k in z.files:
            w = wmap.get(k, 0.0)
            if w <= 0:
                continue
            a = z[k].astype("float32").ravel(); a[~np.isfinite(a)] = 1e4
            b = np.clip(np.searchsorted(edges, a, side="right") - 1, 0, len(edges) - 1)
            hist += np.bincount(b * (H * W) + np.arange(H * W), minlength=hist.size).reshape(hist.shape) * w; tot += w
    cum = np.cumsum(hist, axis=0) / tot
    med_bin = np.argmax(cum >= 0.5, axis=0)
    med = edges[med_bin] + 0.5; med[med_bin >= len(edges) - 1] = np.nan
    med = med.reshape(H, W)
    xs = cfg["grid"]["minx"] + (np.arange(W) + 0.5) * RES; ys = cfg["grid"]["maxy"] - (np.arange(H) + 0.5) * RES
    fig0, ax0 = plt.subplots(); cs = ax0.contour(xs, ys, np.where(np.isfinite(med), med, 999), levels=LEVELS); plt.close(fig0)
    tr = Transformer.from_crs(3400, 4326, always_xy=True)
    from shapely.geometry import LineString, mapping
    feats = []
    for lev, segs in zip(cs.levels, cs.allsegs):
        for sg in segs:
            if len(sg) < 2:
                continue
            lo, la = tr.transform(sg[:, 0], sg[:, 1])
            feats.append({"type": "Feature", "properties": {"minutes_after_1500": float(lev), "clock_mdt": C.clock(lev)},
                          "geometry": mapping(LineString(np.round(np.c_[lo, la], 6)))})
    json.dump({"type": "FeatureCollection", "name": "Rubicon Ridge v5 weighted-median fire arrival (scar-limited runs, default weights)",
               "note": "Per 40 m cell, the weighted median (v5 fit weights) of the arrival time over the 1,440 weighted scar-limited runs. "
                       "Model output, not an observation. Times MDT, Oct. 17, 2026.", "features": feats},
              open(ISO, "w"))
    print("isochrones:", len(feats), "lines")
iso = gpd.read_file(ISO) if os.path.exists(ISO) else None

# ---------------------------------------------------------------- figure 1: tracks map (regional / 40 km / zoom)
cmap = plt.get_cmap("viridis")
tnum = lambda t: (pd.to_datetime(t) - T0).dt.total_seconds() / 3600 if hasattr(t, "dt") else (pd.Timestamp(t) - T0).total_seconds() / 3600
norm = matplotlib.colors.Normalize(0, (T1 - T0).total_seconds() / 3600)
fig, axs = plt.subplots(1, 3, figsize=(17, 6.6), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
asp = 1 / np.cos(np.radians(52.8))
seg = ctx[(ctx.mdt >= T0) & (ctx.mdt <= T1)].sort_values(["reg", "mdt"])
for ax, ext, title in [(axs[0], (-115.0, -106.0, 50.0, 57.0), "Whole flights, 3 p.m. Oct. 17 - 4:30 a.m. Oct. 18"),
                       (axs[1], (-111.55, -110.25, 52.43, 53.17), "Supplied file: within 40 km of 52.80 N 110.90 W, at or below 12,000 ft"),
                       (axs[2], (-111.08, -110.74, 52.70, 52.88), "Near the burn, with v5 modelled fire arrival")]:
    if ax is axs[1]:
        d = S.copy(); d["h"] = (d.mdt - T0).dt.total_seconds() / 3600
        ax.scatter(d.lon, d.lat, c=cmap(norm(d.h)), s=3, lw=0, zorder=3)
        th = np.linspace(0, 2 * np.pi, 361)
        ax.plot(-110.90 + 40 / 111.32 / np.cos(np.radians(52.8)) * np.cos(th), 52.80 + 40 / 110.57 * np.sin(th), color=MUTED, lw=0.7, ls=":")
        ax.text(-110.90, 52.80 + 40 / 110.57 + 0.005, "40 km from the file's 'base' point", ha="center", fontsize=7, color=MUTED)
    else:
        for reg, g in seg.groupby("reg"):
            g = g[(g.lon > ext[0] - 1) & (g.lon < ext[1] + 1)]
            if len(g) < 2:
                continue
            h = ((g.mdt - T0).dt.total_seconds() / 3600).values
            brk = np.r_[False, np.diff(g.mdt.values).astype("timedelta64[s]").astype(float) > 300]
            for i in range(1, len(g)):
                if not brk[i]:
                    ax.plot(g.lon.values[i - 1:i + 1], g.lat.values[i - 1:i + 1], color=cmap(norm(h[i])), lw=1.0 if ax is axs[0] else 1.4, zorder=3)
    burn.boundary.plot(ax=ax, color=FIRE_RED, lw=1.2, zorder=4)
    if ax is axs[2] and iso is not None:
        for lev, g in iso.groupby("minutes_after_1500"):
            g.plot(ax=ax, color=FIRE_RED, lw=0.6, alpha=0.75, ls="--", zorder=4)
            c0 = g.geometry.iloc[int(np.argmax(g.length))].interpolate(0.5, normalized=True)
            ax.text(c0.x, c0.y, C.clock(lev), fontsize=6.5, color=FIRE_RED, ha="center", va="center",
                    bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.8), zorder=5)
        for p in ("Range 2 (RSO)", "Range 15 (RSO)", "Road head"):
            lo, la = cfg["reference_points"][p]
            ax.plot(lo, la, "o", ms=3.5, color=INK, zorder=6)
        ax.annotate("Range 2 / Range 15 / road head", tuple(cfg["reference_points"]["Range 2 (RSO)"]), xytext=(-110.99, 52.80),
                    fontsize=6.5, color=INK, ha="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    marks = [("fire centre (v5)", C.FIRE[1], C.FIRE[0], "*", FIRE_RED, 11)] + \
            [(k.split(" (")[0].replace(" Intl", ""), C.PLACES[k][1], C.PLACES[k][0], "s" if k.startswith("C") else "^", INK, 5)
             for k in (["CYWV Wainwright Airport", "CFB Wainwright airfield (CZWN stn)", "Wainwright town"] if ax is not axs[0] else
                       ["CYEG Edmonton Intl", "CYYC Calgary Intl", "CYXE Saskatoon", "CYLL Lloydminster", "CYLB Lac La Biche",
                        "CYMM Fort McMurray", "CYXH Medicine Hat", "CYWV Wainwright Airport"])]
    for lab, lo, la, mk, col, ms in marks:
        if ext[0] < lo < ext[1] and ext[2] < la < ext[3]:
            ax.plot(lo, la, mk, color=col, ms=ms, mec="white", mew=0.5, zorder=7)
            ax.text(lo + (ext[1] - ext[0]) * 0.012, la + (ext[3] - ext[2]) * 0.008, lab.replace("CFB Wainwright airfield", "CFB airfield (CZWN)"),
                    fontsize=7, color=col, zorder=7, bbox=dict(fc="white", ec="none", pad=0.3, alpha=0.7))
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect(asp)
    ax.set_title(title, fontsize=8.5, loc="left", color=INK)
    ax.tick_params(labelsize=7)
    if ax is axs[0]:
        ax.axvline(-110.005, color=MUTED, lw=0.6); ax.text(-109.95, 50.15, "ALTA. | SASK.", fontsize=6.5, color=MUTED)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
cb = fig.colorbar(sm, ax=axs, orientation="horizontal", fraction=0.04, pad=0.08, aspect=60)
ticks = np.arange(0, 14, 1); cb.set_ticks(ticks); cb.set_ticklabels([(T0 + pd.Timedelta(hours=float(h))).strftime("%H:%M") for h in ticks])
cb.set_label("time of position (MDT); front reached the fire 15:28; first inbound air-ambulance take-off 16:55-16:58", fontsize=8)
axs[2].legend(handles=[Line2D([], [], color=FIRE_RED, lw=1.2, label="burn outline (Oct. 20 vectors)"),
                       Line2D([], [], color=FIRE_RED, lw=0.6, ls="--", label="v5 weighted-median arrival (model)"),
                       Line2D([], [], marker="*", color=FIRE_RED, lw=0, label="fire centre (v5 front-plane origin)")],
              loc="lower left", fontsize=6.5, frameon=True, framealpha=0.85)
fig.suptitle("Rubicon Cloud (v5.1): aircraft heard near CFB Wainwright after the fire, coloured by time", x=0.01, ha="left",
             color=INK, fontsize=12)
fig.text(0.01, 0.005, FOOT, fontsize=6.5, color=MUTED, va="bottom")
fig.subplots_adjust(left=0.03, right=0.99, top=0.9, bottom=0.25, wspace=0.12)
fig.savefig(os.path.join(OUT, "fig_tracks_map.png"), dpi=140); plt.close(fig)

# ---------------------------------------------------------------- figure 2: time-altitude and time-distance per aircraft
regs = AS.sort_values("first_mdt").registration.tolist()
fig, axs = plt.subplots(len(regs), 2, figsize=(15, 1.55 * len(regs) + 1.6), sharex=True)
for i, reg in enumerate(regs):
    g = ctx[(ctx.reg == reg) & (ctx.mdt >= T0) & (ctx.mdt <= T1)].sort_values("mdt")
    s = S[S.registration == reg]
    col = OPS[opcode(cs_of.get(reg, ""))][0]
    for j, (ax, ycol, sy, lim) in enumerate([(axs[i, 0], "alt_ft", "alt_ft", 12000), (axs[i, 1], "km_fire", None, None)]):
        y = g[ycol] if ycol in g else None
        gg = g.copy(); gg.loc[gg.mdt.diff().dt.total_seconds() > 300, ycol] = np.nan     # gaps
        ax.plot(gg.mdt, gg[ycol], color=GRID, lw=0.9, zorder=2)
        if j == 0:
            ax.scatter(s.mdt, s.alt_ft, s=2.5, color=col, zorder=3, lw=0)
            ax.axhline(12000, color=MUTED, lw=0.5, ls=":"); ax.set_ylim(0, 26000)
            ax.set_yticks([0, 10000, 20000]); ax.set_yticklabels(["0", "10k", "20k"], fontsize=7)
        else:
            kf = C.hav(s.lat.values, s.lon.values, *C.FIRE)
            ax.scatter(s.mdt, kf, s=2.5, color=col, zorder=3, lw=0)
            ax.set_yscale("log"); ax.set_ylim(2, 600); ax.set_yticks([3, 10, 30, 100, 300]); ax.set_yticklabels(["3", "10", "30", "100", "300"], fontsize=7)
        ax.axvline(FRONT, color=FIRE_RED, lw=0.9)
        for lab, w in WIN.items():
            ax.axvspan(hm(w["p10"]), hm(w["p90"]), color=FIRE_RED, alpha=0.10, lw=0)
        ax.axvline(pd.Timestamp("2026-10-17 16:00"), color=INK, lw=0.6, ls="--")
        ax.grid(axis="y", color=GRID, lw=0.4); ax.spines[["top", "right"]].set_visible(False)
        for v in V[V.reg == reg].itertuples():
            if pd.notna(v.est_on_ground_mdt):
                ax.axvline(v.est_on_ground_mdt, color=col, lw=0.6, alpha=0.6, ymax=0.25 if j == 0 else 1)
    cs = cs_of.get(reg, ""); tp = AS.loc[AS.registration == reg, "type"].iloc[0]
    axs[i, 0].set_ylabel(f"{reg}\n{tp} {cs}", fontsize=7.5, rotation=0, ha="right", va="center", color=col)
axs[0, 0].set_title("Barometric altitude, ft (grey: whole track from the archive; dots: positions in the supplied file)", fontsize=8.5, loc="left", color=INK)
axs[0, 1].set_title("Distance from the fire centre, km (log scale)", fontsize=8.5, loc="left", color=INK)
for ax in axs[-1]:
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=1)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlim(T0, T1); ax.tick_params(axis="x", labelsize=7, rotation=0)
handles = [Line2D([], [], color=FIRE_RED, lw=0.9, label="front at the fire 15:28 (v5)"),
           Patch(color=FIRE_RED, alpha=0.25, label=f"v5 arrival 10-90%: Range 2 {WIN['Range 2']['p10']}-{WIN['Range 2']['p90']}, "
                                                   f"Range 15 {WIN['Range 15']['p10']}-{WIN['Range 15']['p90']}, road head {WIN['road head']['p10']}-{WIN['road head']['p90']}"),
           Line2D([], [], color=INK, lw=0.6, ls="--", label="16:00 (DND: incident 'just before 4:00')"),
           Line2D([], [], color=MUTED, lw=0.6, alpha=0.8, label="estimated landing / take-off, Wainwright area (coloured ticks)")] + \
          [Line2D([], [], marker="o", lw=0, color=c, label=l) for c, l in OPS.values()]
fig.legend(handles=handles, loc="upper left", ncol=3, fontsize=7, frameon=False, bbox_to_anchor=(0.01, 0.965))
fig.suptitle("Rubicon Cloud (v5.1): each aircraft's altitude and distance from the fire against the fire timeline", x=0.01, ha="left", color=INK, fontsize=12)
fig.text(0.01, 0.003, FOOT, fontsize=6.5, color=MUTED, va="bottom")
fig.subplots_adjust(left=0.1, right=0.99, top=0.89, bottom=0.07, hspace=0.25, wspace=0.08)
fig.savefig(os.path.join(OUT, "fig_time_alt_dist.png"), dpi=130); plt.close(fig)

# ---------------------------------------------------------------- figure 3: response timeline (Gantt) with a zoom on 15:00-18:30
air = json.load(open(os.path.join(OUT, "air_timing_check.json")))
gap = air["airspace_null_check"]["gap_below_12000ft_mdt"]
fig, axs = plt.subplots(2, 1, figsize=(15, 9.5), gridspec_kw=dict(height_ratios=[1, 1.25]))
order = regs[::-1]
for ax, (a, b) in zip(axs, [(T0, pd.Timestamp("2026-10-17 18:45")), (T0, T1)]):
    for yi, reg in enumerate(order):
        col = OPS[opcode(cs_of.get(reg, ""))][0]
        for f in FL[FL.reg == reg].itertuples():
            s0 = f.departure_mdt if pd.notna(f.departure_mdt) else f.first_heard_mdt
            s1 = f.arrival_mdt if pd.notna(f.arrival_mdt) else f.last_heard_mdt
            if s1 < a or s0 > b:
                continue
            ax.plot([f.first_heard_mdt, f.last_heard_mdt], [yi, yi], color=col, lw=5, solid_capstyle="butt", zorder=3)
            ax.plot([s0, f.first_heard_mdt], [yi, yi], color=col, lw=5, alpha=0.35, solid_capstyle="butt", zorder=3)
            ax.plot([f.last_heard_mdt, s1], [yi, yi], color=col, lw=5, alpha=0.35, solid_capstyle="butt", zorder=3)
            if ax is axs[0]:
                lab = f"{str(f.origin).split(' ')[0]}" + "→" + f"{str(f.destination).split(' ')[0]}"
                ax.text(f.first_heard_mdt, yi + 0.3, lab, fontsize=6.5, color=col, va="bottom", clip_on=True)
        for v in V[V.reg == reg].itertuples():
            if pd.notna(v.est_on_ground_mdt) and a <= v.est_on_ground_mdt <= b:
                ax.plot(v.est_on_ground_mdt, yi, marker="v" if v.event == "arrival" else "^", color=INK, ms=5, zorder=5)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{r} {cs_of.get(r, '')} ({AS.loc[AS.registration == r, 'type'].iloc[0]})" for r in order], fontsize=7.5)
    ax.axvline(FRONT, color=FIRE_RED, lw=1)
    for lab, w in WIN.items():
        ax.axvspan(hm(w["p10"]), hm(w["p90"]), color=FIRE_RED, alpha=0.08, lw=0)
    ax.axvline(pd.Timestamp("2026-10-17 16:00"), color=INK, lw=0.7, ls="--")
    ax.axvspan(hm(gap[0]), hm(gap[1]), ymin=0, ymax=1, color=GRID, alpha=0.35, lw=0, zorder=0)
    ax.set_xlim(a, b); ax.grid(axis="x", color=GRID, lw=0.4); ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 30]) if ax is axs[0] else mdates.HourLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M")); ax.tick_params(axis="x", labelsize=7.5)
s, e = FR.iloc[0], FR[FR.callsigns.fillna("").str.startswith("STAR")].iloc[0]
ax = axs[0]
ytop = len(order) - 0.4
ax.annotate(f"front at the fire\n15:28 (v5)", (FRONT, ytop), xytext=(4, 0), textcoords="offset points", fontsize=7, color=FIRE_RED, va="top")
ax.annotate("v5 arrival windows\n(10-90%) at Range 2,\nRange 15, road head", (hm(WIN["Range 2"]["p10"]), 1.0), xytext=(4, 0),
            textcoords="offset points", fontsize=7, color=FIRE_RED)
ax.annotate("DND: 'just before 4:00'", (pd.Timestamp("2026-10-17 16:00"), ytop - 1.4), xytext=(4, 0), textcoords="offset points", fontsize=7, color=INK)
ax.annotate(f"no aircraft below 12,000 ft heard within 80 km of the fire\n{gap[0][:5]}-{gap[1][:5]} ({gap[2]:.0f} min)",
            (0.31, 0.40), xycoords="axes fraction", fontsize=7, color=MUTED, ha="left")
ax.annotate(f"first inbound take-offs:\n{s.reg} {pd.Timestamp(s.departure_or_first_heard_mdt):%H:%M} (estimated, field below coverage)\n"
            f"{e.reg} {e.callsigns} {pd.Timestamp(e.departure_or_first_heard_mdt):%H:%M} (heard, Edmonton Intl)",
            (0.31, 0.10), xycoords="axes fraction", fontsize=7, color=INK)
axs[0].set_title("3:00-6:45 p.m.: solid = heard; faint = estimated before first / after last position; ▼ landing, ▲ take-off in the "
                 "Wainwright area (estimated)", fontsize=8.5, loc="left", color=INK)
axs[1].set_title("The whole evening and night: repeated shuttles between the Wainwright area and Edmonton, Calgary and Saskatoon",
                 fontsize=8.5, loc="left", color=INK)
axs[1].legend(handles=[Patch(color=c, label=l) for c, l in OPS.values()], loc="lower right", fontsize=7, frameon=True, framealpha=0.9)
fig.suptitle("Rubicon Cloud (v5.1): the air response in time, against the front and the v5 fire-arrival windows", x=0.01, ha="left",
             color=INK, fontsize=12)
fig.text(0.01, 0.003, FOOT, fontsize=6.5, color=MUTED, va="bottom")
fig.subplots_adjust(left=0.15, right=0.99, top=0.92, bottom=0.09, hspace=0.25)
fig.savefig(os.path.join(OUT, "fig_response_timeline.png"), dpi=140); plt.close(fig)

# ---------------------------------------------------------------- optional animation
if args.gif:
    from PIL import Image
    import io
    ext = (-111.75, -110.05, 52.35, 53.25)
    d = ctx[(ctx.lon > ext[0]) & (ctx.lon < ext[1]) & (ctx.lat > ext[2]) & (ctx.lat < ext[3])].sort_values("mdt")
    frames = []
    times = pd.date_range(T0, pd.Timestamp("2026-10-18 04:00"), freq="10min")
    for t in times:
        fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=80)
        burn.boundary.plot(ax=ax, color=FIRE_RED, lw=1.0)
        trail = d[(d.mdt > t - pd.Timedelta(minutes=30)) & (d.mdt <= t)]
        for reg, g in trail.groupby("reg"):
            col = OPS[opcode(cs_of.get(reg, ""))][0]
            ax.plot(g.lon, g.lat, color=col, lw=1.2, alpha=0.8)
            ax.plot(g.lon.iloc[-1], g.lat.iloc[-1], "o", color=col, ms=4)
            ax.text(g.lon.iloc[-1] + 0.02, g.lat.iloc[-1], f"{reg} {g.alt_baro.iloc[-1]} ft", fontsize=6, color=col)
        for k in ["CYWV Wainwright Airport", "CFB Wainwright airfield (CZWN stn)", "Wainwright town"]:
            ax.plot(C.PLACES[k][1], C.PLACES[k][0], "s", color=INK, ms=3)
        ax.plot(C.FIRE[1], C.FIRE[0], "*", color=FIRE_RED, ms=9)
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect(asp); ax.tick_params(labelsize=6)
        st = "front at the fire 15:28" if t < pd.Timestamp("2026-10-17 16:00") else ""
        ax.set_title(f"{t:%b %d %H:%M} MDT   (last 30 min of tracks)   {st}", fontsize=8, loc="left", color=INK)
        fig.text(0.01, 0.01, "adsb.lol (ODbL); burn outline: Prime Rogue vectors; Prime Rogue Inc.", fontsize=5.5, color=MUTED)
        buf = io.BytesIO(); fig.savefig(buf, format="png"); plt.close(fig); buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=64))
    frames[0].save(os.path.join(OUT, "anim_tracks_1500-0400.gif"), save_all=True, append_images=frames[1:], duration=250, loop=0, optimize=True)
    print("gif bytes", os.path.getsize(os.path.join(OUT, "anim_tracks_1500-0400.gif")))
print("figures done")
