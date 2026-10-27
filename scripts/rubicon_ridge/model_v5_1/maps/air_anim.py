"""Map and animate the low-altitude aircraft positions near CFB Wainwright, Oct. 17-18, 2026 (MDT).
Usage: python3 air_anim.py <positions.csv> <out_dir>"""
import sys, os
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.lines import Line2D
from PIL import Image

CSV, OUT = sys.argv[1], sys.argv[2]
B = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../.."))
A = os.environ.get("CFB_ARCHIVE_1", os.path.join(B, "../cfb-wainwright-fire-archive"))
LOGO = os.environ.get("PR_LOGO", "")   # optional Prime Rogue badge PNG

d = pd.read_csv(CSV, parse_dates=["mdt"]).sort_values("mdt")
W, E, S, N = -111.52, -110.27, 52.50, 53.19
LAT0 = 52.85; ASP = 1 / np.cos(np.radians(LAT0))

# ---- base layers
osm = gpd.read_file(f"{A}/data/infrastructure/osm_infrastructure_ODbL.geojson").cx[W:E, S:N]
water = osm[osm.kind == "water"]
roads = osm[(osm.kind == "road") & osm.highway.isin(["motorway", "trunk", "primary", "secondary"])]
base = osm[osm.name == "CFB Wainwright"]
pda = osm[osm.name.fillna("").str.startswith("Permanent Danger Area")]
burn = gpd.read_file(f"{B}/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson").to_crs(4326)
PLACES = {"Wainwright": (-110.8601, 52.8334), "Edgerton": (-110.448, 52.749), "Irma": (-111.234, 52.913),
          "Fabyan": (-110.973, 52.867)}
FIRE = (-110.9052, 52.7645); CZWN = (-111.096, 52.829)
FACS = [("Wainwright Airport (CYWV)\nrunway 13/31, 3,000 ft paved", (-110.8587, 52.7967), "s", 40, -38),
        ("Field 21 airfield (CZWN)", (-111.096, 52.829), "s", -70, -48),
        ("Camp Wainwright helipad (CFF7)", (-110.9045, 52.8280), "H", -40, 48),
        ("Health Centre heliport (CWH2)", (-110.8647, 52.8430), "H", 30, 34)]
FACS = [("Wainwright Airport (CYWV)\nrunway 13/31, 3,000 ft paved", (-110.8587, 52.7967), "s", 40, -38),
        ("Field 21 airfield (CZWN)", (-111.096, 52.829), "s", -70, -48),
        ("Camp Wainwright helipad (CFF7)", (-110.9045, 52.8280), "H", -40, 48),
        ("Health Centre heliport (CWH2)", (-110.8647, 52.8430), "H", 30, 34)]

# ---- aircraft colours
regs = (d.groupby("registration").mdt.min().sort_values().index.tolist())
desc = d.groupby("registration").description.first()
PAL = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6",
       "#bfef45", "#fabed4", "#ffffff"]
COL = {r: PAL[i % len(PAL)] for i, r in enumerate(regs)}
def short(s):
    return (s.replace("AIRBUS HELICOPTERS ", "Airbus ").replace("EUROCOPTER ", "Eurocopter ")
             .replace("BEECHCRAFT ", "Beech ").replace("BEECH ", "Beech ").replace("CESSNA ", "Cessna "))
LBL = {r: f"{r}  {short(desc[r])}" for r in regs}

BG, LAND, INK = "#0f1418", "#1b2329", "#e9edf0"
def draw_base(ax):
    ax.set_facecolor(LAND)
    base.plot(ax=ax, facecolor="#2a3a2a", edgecolor="#6f8f6f", lw=0.8, zorder=1)
    pda.plot(ax=ax, facecolor="none", edgecolor="#8a7a4a", lw=0.5, ls="--", zorder=1)
    water.plot(ax=ax, facecolor="#24465e", edgecolor="none", zorder=2)
    roads.plot(ax=ax, color="#59636e", lw=0.6, zorder=2)
    burn.plot(ax=ax, facecolor="#7a1f12", edgecolor="#ff6a3d", lw=0.8, alpha=0.9, zorder=3)
    ax.plot(*FIRE, marker="*", color="#ff6a3d", ms=11, mec="white", mew=0.6, zorder=6)
    ax.text(FIRE[0] + 0.015, FIRE[1] - 0.018, "Oct. 17 burn", color="#ff9d7a", fontsize=8, zorder=6)
    for n, (x, y) in PLACES.items():
        ax.plot(x, y, "o", color=INK, ms=3, zorder=6); ax.text(x + 0.01, y + 0.005, n, color=INK, fontsize=8, zorder=6)
    ax.text(-110.99, 52.69, "CFB Wainwright", color="#8fb08f", fontsize=8, style="italic", zorder=6)
    # ---- airstrips / heliports and distance rings from the fire centre
    kxk = 111.32 * np.cos(np.radians(LAT0)); kyk = 110.95
    for rk in (5, 10, 20, 30, 40):
        th = np.linspace(0, 2*np.pi, 361)
        ax.plot(FIRE[0] + rk/kxk*np.sin(th), FIRE[1] + rk/kyk*np.cos(th), color="#ff9d7a", lw=0.6, ls=(0, (4, 4)), alpha=0.55, zorder=3)
        ax.text(FIRE[0] + rk/kxk*np.sin(np.radians(225)), FIRE[1] + rk/kyk*np.cos(np.radians(225)), f"{rk} km",
                color="#ff9d7a", fontsize=7, ha="center", va="center", zorder=6, clip_on=True,
                bbox=dict(facecolor=LAND, edgecolor="none", pad=0.5))
    # CYWV runway 13/31: 3,000 ft (914 m), true heading 127 deg, centre 52.7967 N 110.8587 W
    hl = 0.457; hd = np.radians(127); cx, cy = -110.8587, 52.7967
    ax.plot([cx - hl/kxk*np.sin(hd), cx + hl/kxk*np.sin(hd)], [cy - hl/kyk*np.cos(hd), cy + hl/kyk*np.cos(hd)],
            color="white", lw=4, solid_capstyle="butt", zorder=7)
    ax.plot([cx - hl/kxk*np.sin(hd), cx + hl/kxk*np.sin(hd)], [cy - hl/kyk*np.cos(hd), cy + hl/kyk*np.cos(hd)],
            color="black", lw=1, ls=(0, (2, 2)), zorder=7)
    for name, (x, y), m, dx, dy in FACS:
        dk = np.hypot((x - FIRE[0])*kxk, (y - FIRE[1])*kyk)
        ax.plot(x, y, m, color="white", mec="black", ms=7, zorder=7)
        ax.annotate(f"{name}\n{dk:.1f} km from fire", (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=7, color=INK,
                    arrowprops=dict(arrowstyle="-", color="#8b949e", lw=0.6), zorder=7,
                    bbox=dict(facecolor="#0f1418", edgecolor="#30363d", alpha=0.85, pad=1.5))
    # ---- airstrips / heliports and distance rings from the fire centre
    kxk = 111.32 * np.cos(np.radians(LAT0)); kyk = 110.95
    for rk in (5, 10, 20, 30, 40):
        th = np.linspace(0, 2*np.pi, 361)
        ax.plot(FIRE[0] + rk/kxk*np.sin(th), FIRE[1] + rk/kyk*np.cos(th), color="#ff9d7a", lw=0.6, ls=(0, (4, 4)), alpha=0.55, zorder=3)
        ax.text(FIRE[0] + rk/kxk*np.sin(np.radians(225)), FIRE[1] + rk/kyk*np.cos(np.radians(225)), f"{rk} km",
                color="#ff9d7a", fontsize=7, ha="center", va="center", zorder=6, clip_on=True,
                bbox=dict(facecolor=LAND, edgecolor="none", pad=0.5))
    # CYWV runway 13/31: 3,000 ft (914 m), true heading 127 deg, centre 52.7967 N 110.8587 W
    hl = 0.457; hd = np.radians(127); cx, cy = -110.8587, 52.7967
    xs = [cx - hl/kxk*np.sin(hd), cx + hl/kxk*np.sin(hd)]; ys = [cy - hl/kyk*np.cos(hd), cy + hl/kyk*np.cos(hd)]
    ax.plot(xs, ys, color="white", lw=4, solid_capstyle="butt", zorder=7)
    ax.plot(xs, ys, color="black", lw=1, ls=(0, (2, 2)), zorder=7)
    for name, (x, y), m, dx, dy in FACS:
        dk = np.hypot((x - FIRE[0])*kxk, (y - FIRE[1])*kyk)
        ax.plot(x, y, m, color="white", mec="black", ms=7, zorder=7)
        ax.annotate(f"{name}\n{dk:.1f} km from fire", (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=7, color=INK,
                    arrowprops=dict(arrowstyle="-", color="#8b949e", lw=0.6), zorder=7,
                    bbox=dict(facecolor="#0f1418", edgecolor="#30363d", alpha=0.85, pad=1.5))
    # 10 km scale bar
    kx = 111.32 * np.cos(np.radians(LAT0)); x0, y0 = W + 0.04, S + 0.03
    ax.plot([x0, x0 + 10 / kx], [y0, y0], color=INK, lw=2.5, zorder=7); ax.text(x0, y0 + 0.012, "10 km", color=INK, fontsize=8, zorder=7)
    ax.set_xlim(W, E); ax.set_ylim(S, N); ax.set_aspect(ASP)
    ax.tick_params(colors="#8b949e", labelsize=7)
    for s in ax.spines.values(): s.set_color("#30363d")

def legend(ax, loc="upper right"):
    h = [Line2D([], [], color=COL[r], lw=2.5, marker="o", ms=5, label=LBL[r]) for r in regs]
    h += [Line2D([], [], ls="", marker="*", color="#ff6a3d", ms=10, label="Oct. 17 burn (fire centre)"),
          Line2D([], [], ls="", marker="s", color="white", mec="black", ms=7, label="airfield (fixed wing)"),
          Line2D([], [], ls="", marker="H", color="white", mec="black", ms=7, label="helipad / heliport"),
          Line2D([], [], color="#ff9d7a", ls=(0, (4, 4)), label="distance from fire centre (km)")]
    lg = ax.legend(handles=h, loc=loc, fontsize=7, frameon=True, facecolor="#0f1418", edgecolor="#30363d",
                   labelcolor=INK, title="Aircraft (registration, type)", title_fontsize=7.5)
    lg.get_title().set_color(INK)

def logo(fig, ax):
    if os.path.exists(LOGO):
        ia = ax.inset_axes([0.86, 0.015, 0.13, 0.13], zorder=10)
        ia.imshow(Image.open(LOGO)); ia.axis("off")

CREDIT = ("Prime Rogue Inc. | Aircraft: ADS-B positions below about 12,000 ft within 40 km of the base (no military transponders in this set) | "
          "Burn outline: Prime Rogue, from Oct. 20 imagery | Map data © OpenStreetMap contributors (ODbL)")

# ---- static map
fig, ax = plt.subplots(figsize=(12, 10.5), dpi=130, facecolor=BG)
draw_base(ax)
for r in regs:
    s = d[d.registration == r]
    ax.plot(s.lon, s.lat, color=COL[r], lw=1.2, alpha=0.9, zorder=4)
    ax.plot(s.lon.iloc[0], s.lat.iloc[0], "o", color=COL[r], ms=5, mec="black", mew=0.5, zorder=5)
legend(ax)
ax.set_title("Aircraft near CFB Wainwright after the fire, Oct. 17, 3:31 p.m. to Oct. 18, 3:58 a.m. (MDT)\n"
             "Full tracks; dot marks first position", color=INK, fontsize=11, loc="left")
fig.text(0.01, 0.01, CREDIT, color="#8b949e", fontsize=6.5)
logo(fig, ax)
fig.savefig(f"{OUT}/aircraft_tracks_map.png", facecolor=BG, bbox_inches="tight")
plt.close(fig)

# ---- animation
STEP = pd.Timedelta(minutes=2); TRAIL = pd.Timedelta(minutes=20); GAP = pd.Timedelta(minutes=10)
times = pd.date_range(d.mdt.min().floor("2min"), d.mdt.max().ceil("2min"), freq=STEP)
fig, ax = plt.subplots(figsize=(10, 9), dpi=100, facecolor=BG)
draw_base(ax); legend(ax)
lines = {r: ax.plot([], [], color=COL[r], lw=2.0, alpha=0.95, zorder=4)[0] for r in regs}
dots = {r: ax.plot([], [], "o", color=COL[r], ms=7, mec="black", mew=0.7, zorder=5)[0] for r in regs}
tags = {r: ax.text(0, 0, "", color=COL[r], fontsize=7, weight="bold", zorder=6) for r in regs}
clock = ax.text(0.015, 0.975, "", transform=ax.transAxes, color=INK, fontsize=15, weight="bold", va="top",
                bbox=dict(facecolor="#0f1418", edgecolor="#30363d", pad=5), zorder=8)
ax.set_title("Aircraft near CFB Wainwright after the Oct. 17, 2026 fire (MDT)", color=INK, fontsize=11, loc="left")
fig.text(0.01, 0.008, CREDIT, color="#8b949e", fontsize=5.5)
logo(fig, ax)
G = {r: d[d.registration == r] for r in regs}
def frame(i):
    t = times[i]
    for r in regs:
        s = G[r]; tr = s[(s.mdt <= t) & (s.mdt > t - TRAIL)]
        lines[r].set_data(tr.lon, tr.lat)
        last = s[s.mdt <= t].tail(1)
        if len(last) and t - last.mdt.iloc[0] <= GAP:
            x, y = last.lon.iloc[0], last.lat.iloc[0]
            dots[r].set_data([x], [y]); tags[r].set_position((x + 0.012, y + 0.006))
            tags[r].set_text(f"{r}\n{int(last.alt_ft.iloc[0]):,} ft")
        else:
            dots[r].set_data([], []); tags[r].set_text("")
    day = "Oct. 17" if t.day == 17 else "Oct. 18"
    clock.set_text(f"{day}  {t.strftime('%-I:%M %p').lower().replace('am','a.m.').replace('pm','p.m.')}")
    return list(lines.values()) + list(dots.values()) + list(tags.values()) + [clock]
anim = FuncAnimation(fig, frame, frames=len(times), blit=False)
anim.save(f"{OUT}/aircraft_tracks_animation.mp4", writer=FFMpegWriter(fps=12, bitrate=2500), savefig_kwargs={"facecolor": BG})
print("frames", len(times))
