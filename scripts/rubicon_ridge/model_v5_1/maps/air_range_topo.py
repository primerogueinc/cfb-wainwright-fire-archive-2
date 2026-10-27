"""Aircraft tracks on a range map with topography (MRDEM-30 hillshade), CFB Wainwright, Oct. 17-18, 2026.
Usage: python3 air_range_topo.py <positions.csv> <out_png>"""
import sys, os
import numpy as np, pandas as pd, geopandas as gpd, rasterio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D
from shapely.geometry import box
from PIL import Image

CSV, OUTPNG = sys.argv[1], sys.argv[2]
B = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../..")); A = os.environ.get("CFB_ARCHIVE_1", os.path.join(B, "../cfb-wainwright-fire-archive"))
LOGO = os.environ.get("PR_LOGO", "")   # optional Prime Rogue badge PNG
CRS = "EPSG:3400"

dem = rasterio.open(f"{A}/data/dem/mrdem30_dtm_3400_10m.tif")
z = dem.read(1, masked=True).astype(float).filled(np.nan)
x0, y0, x1, y1 = dem.bounds.left, dem.bounds.bottom, dem.bounds.right, dem.bounds.top
ls = LightSource(azdeg=315, altdeg=40)
valid = np.isfinite(z) & (z > 0)
z[~valid] = np.nan
zz = np.where(valid, z, np.nanmedian(z))
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("t", plt.get_cmap("terrain")(np.linspace(0.25, 0.95, 256)))
rgb = ls.shade(zz, cmap=cmap, blend_mode="soft", vert_exag=6, vmin=np.nanpercentile(z, 1), vmax=np.nanpercentile(z, 99.5))
rgb = np.dstack([rgb[..., :3], np.where(valid, 1.0, 0.0)])

def to3400(g): return g.to_crs(CRS)
ext = box(x0, y0, x1, y1)
osm = to3400(gpd.read_file(f"{A}/data/infrastructure/osm_infrastructure_ODbL.geojson"))
osm = osm[osm.intersects(ext)]
water = osm[osm.kind == "water"]
roads = osm[(osm.kind == "road")]
mil = osm[osm.kind == "military"]
cfb = mil[mil.name == "CFB Wainwright"]
pda = mil[mil.name.fillna("").str.contains("Danger Area")]
rng = mil[~mil.name.isin(["CFB Wainwright"]) & ~mil.name.fillna("").str.contains("Danger Area")]
burn = to3400(gpd.read_file(f"{B}/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson"))
rso = pd.read_csv(f"{A}/data/reference/range_standing_orders_2024/rso_points.csv")
rso = rso[~rso.name.str.contains("barrier")]
rsog = to3400(gpd.GeoDataFrame(rso, geometry=gpd.points_from_xy(rso.lon, rso.lat), crs=4326))

FAC = {"Wainwright Airport (CYWV)": (-110.8587, 52.7967), "Camp Wainwright helipad (CFF7)": (-110.9045, 52.8280),
       "Wainwright Health Centre heliport (CWH2)": (-110.8647, 52.8430), "Field 21 airfield (CZWN)": (-111.096, 52.829)}
facg = to3400(gpd.GeoDataFrame({"name": list(FAC)}, geometry=gpd.points_from_xy(*zip(*FAC.values())), crs=4326))
TOWN = to3400(gpd.GeoSeries(gpd.points_from_xy([-110.8601], [52.8334]), crs=4326)).iloc[0]
FIRE = to3400(gpd.GeoSeries(gpd.points_from_xy([-110.9052], [52.7645]), crs=4326)).iloc[0]

d = pd.read_csv(CSV, parse_dates=["mdt"]).sort_values("mdt")
g = to3400(gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d.lon, d.lat), crs=4326)); d["x"] = g.geometry.x; d["y"] = g.geometry.y
regs = d.groupby("registration").mdt.min().sort_values().index.tolist()
PAL = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#ffffff"]
COL = {r: PAL[i % len(PAL)] for i, r in enumerate(regs)}
desc = d.groupby("registration").description.first()
def short(s):
    return (s.replace("AIRBUS HELICOPTERS ", "Airbus ").replace("EUROCOPTER ", "Eurocopter ")
             .replace("BEECHCRAFT ", "Beech ").replace("BEECH ", "Beech ").replace("CESSNA ", "Cessna "))

fig, ax = plt.subplots(figsize=(13, 13), dpi=130, facecolor="#0f1418")
ax.set_facecolor("#d9d4c7"); ax.imshow(rgb, extent=(x0, x1, y0, y1), origin="upper", zorder=0)
cs = ax.contour(np.linspace(x0, x1, z.shape[1]), np.linspace(y1, y0, z.shape[0]), zz, levels=np.arange(600, 800, 10)) if False else ax.contour(np.linspace(x0, x1, z.shape[1]), np.linspace(y1, y0, z.shape[0]), np.ma.masked_invalid(z), levels=np.arange(600, 800, 10),
                colors="black", linewidths=0.25, alpha=0.35, zorder=1)
water.plot(ax=ax, facecolor="#2f6f9a", edgecolor="none", alpha=0.85, zorder=2)
roads.plot(ax=ax, color="#3a3a3a", lw=0.35, alpha=0.7, zorder=2)
cfb.boundary.plot(ax=ax, color="#1f5c2e", lw=1.6, zorder=3)
pda.boundary.plot(ax=ax, color="#8a5a00", lw=1.0, ls="--", zorder=3)
rng.plot(ax=ax, facecolor="#ffd54a", edgecolor="#7a5c00", lw=0.7, alpha=0.55, zorder=3)
for _, r in rng.iterrows():
    if r["name"]:
        p = r.geometry.representative_point(); ax.text(p.x, p.y, r["name"], fontsize=6, color="#3b2a00", ha="center", zorder=4)
for _, r in pda.iterrows():
    p = r.geometry.representative_point(); ax.text(p.x, p.y, r["name"].replace("Permanent Danger Area", "PDA"), fontsize=8, color="#6b4500", style="italic", zorder=4)
burn.plot(ax=ax, facecolor="#b3261e", edgecolor="#ff6a3d", lw=0.8, alpha=0.75, zorder=4)
ax.plot(FIRE.x, FIRE.y, "*", color="#ff6a3d", ms=13, mec="black", mew=0.6, zorder=8)
for _, r in rsog.iterrows():
    ax.plot(r.geometry.x, r.geometry.y, "D", color="black", mfc="#ffd54a", ms=4.5, zorder=7)
    ax.text(r.geometry.x + 150, r.geometry.y + 120, r["name"], fontsize=6.5, color="black", weight="bold", zorder=7)
for _, r in facg.iterrows():
    ax.plot(r.geometry.x, r.geometry.y, "s", color="white", mec="black", ms=8, zorder=8)
    dk = r.geometry.distance(FIRE) / 1000
    ax.text(r.geometry.x + 250, r.geometry.y - 650, f"{r['name']}\n{dk:.1f} km from fire", fontsize=7.5, color="black", weight="bold",
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.2), zorder=8)
for rk in (2, 5, 10, 15):
    th = np.linspace(0, 2*np.pi, 361); ax.plot(FIRE.x + rk*1000*np.sin(th), FIRE.y + rk*1000*np.cos(th), color="#b3261e", lw=0.9, ls=(0, (5, 4)), alpha=0.8, zorder=4)
    ax.text(FIRE.x + rk*1000*np.sin(np.radians(200)), FIRE.y + rk*1000*np.cos(np.radians(200)), f"{rk} km from fire", fontsize=7.5, color="#7a120c", weight="bold", ha="center",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0.8), zorder=8, clip_on=True)
rw = to3400(gpd.GeoSeries(gpd.points_from_xy([-110.8587], [52.7967]), crs=4326)).iloc[0]; hd = np.radians(127 + 0.0)
ax.plot([rw.x - 457*np.sin(hd), rw.x + 457*np.sin(hd)], [rw.y - 457*np.cos(hd), rw.y + 457*np.cos(hd)], color="white", lw=5, solid_capstyle="butt", zorder=8)
ax.plot([rw.x - 457*np.sin(hd), rw.x + 457*np.sin(hd)], [rw.y - 457*np.cos(hd), rw.y + 457*np.cos(hd)], color="black", lw=1, ls=(0, (2, 2)), zorder=8)
ax.text(TOWN.x + 300, TOWN.y + 250, "Wainwright", fontsize=10, color="black", weight="bold", zorder=8)

# tracks + data-gap events (lost = dropped below coverage; reappear = climbing back into coverage)
for r in regs:
    s = d[d.registration == r].reset_index(drop=True)
    seg = (s.mdt.diff().dt.total_seconds() / 60 > 8).cumsum()
    for _, t in s.groupby(seg):
        ax.plot(t.x, t.y, color=COL[r], lw=2.0, alpha=0.95, zorder=5, solid_capstyle="round")
        ax.plot(t.x, t.y, color="black", lw=3.2, alpha=0.35, zorder=4.9)
        a, b = t.iloc[0], t.iloc[-1]
        for p, m in ((a, "^"), (b, "v")):
            if x0 < p.x < x1 and y0 < p.y < y1:
                ax.plot(p.x, p.y, m, color=COL[r], mec="black", ms=8, zorder=6)
                ax.text(p.x + 200, p.y + 200, f"{p.mdt:%-H:%M}", fontsize=6.5, color="white",
                        bbox=dict(facecolor="black", alpha=0.55, edgecolor="none", pad=0.8), zorder=6)

# scale bar, north arrow
sx, sy = x0 + 1200, y0 + 1200
ax.plot([sx, sx + 5000], [sy, sy], color="black", lw=4, zorder=9); ax.plot([sx, sx + 5000], [sy, sy], color="white", lw=2, zorder=9)
ax.text(sx, sy + 350, "5 km", fontsize=9, weight="bold", zorder=9)
ax.annotate("N", xy=(x1 - 1200, y1 - 900), xytext=(x1 - 1200, y1 - 2400), ha="center", fontsize=12, weight="bold",
            arrowprops=dict(arrowstyle="-|>", color="black", lw=2), zorder=9)
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])

h = [Line2D([], [], color=COL[r], lw=2.5, label=f"{r}  {short(desc[r])}") for r in regs]
h += [Line2D([], [], ls="", marker="v", color="grey", mec="black", ms=8, label="last position before dropping below coverage (time MDT)"),
      Line2D([], [], ls="", marker="^", color="grey", mec="black", ms=8, label="first position climbing back into coverage"),
      Line2D([], [], ls="", marker="s", color="white", mec="black", ms=8, label="airfield / heliport"),
      Line2D([], [], ls="", marker="D", color="black", mfc="#ffd54a", ms=6, label="range (official grid, 2024 Range Standing Orders)"),
      plt.Rectangle((0, 0), 1, 1, facecolor="#ffd54a", alpha=0.55, edgecolor="#7a5c00", label="range area (OpenStreetMap)"),
      Line2D([], [], color="#8a5a00", ls="--", label="Permanent Danger Area"),
      Line2D([], [], color="#1f5c2e", lw=1.6, label="CFB Wainwright boundary"),
      plt.Rectangle((0, 0), 1, 1, facecolor="#b3261e", alpha=0.75, label="Oct. 17 burn"),
      Line2D([], [], ls="", marker="*", color="#ff6a3d", mec="black", ms=11, label="fire centre"),
      Line2D([], [], color="#b3261e", ls=(0, (5, 4)), label="distance from fire centre"),
      Line2D([], [], color="white", lw=4, label="Wainwright Airport runway 13/31 (3,000 ft, paved)")]
lg = ax.legend(handles=h, loc="lower right", fontsize=7, framealpha=0.88, title="Aircraft and map key", title_fontsize=8)
ax.set_title("Air ambulance and helicopter tracks over the CFB Wainwright ranges, Oct. 17, 5:25 p.m. to Oct. 18, 4 a.m. (MDT)\n"
             "Positions exist only above about 4,250 ft (about 2,000 ft above ground): landings and anything flown lower are not visible",
             color="#e9edf0", fontsize=10.5, loc="left")
fig.text(0.01, 0.005, "Prime Rogue Inc. | Aircraft: ADS-B positions (civil transponders only) | Terrain: NRCan MRDEM-30 (Open Government Licence - Canada) | "
         "Ranges/roads/water: © OpenStreetMap contributors (ODbL); range grids: 2024 Range Standing Orders | Burn: Prime Rogue, from Oct. 20 imagery | "
         "Airfields: Transport Canada data via SkyVector/AirportRoutes", color="#8b949e", fontsize=6)
if os.path.exists(LOGO):
    ia = ax.inset_axes([0.005, 0.86, 0.12, 0.12], zorder=10); ia.imshow(Image.open(LOGO)); ia.axis("off")
[t.set_clip_on(True) for t in ax.texts]
fig.savefig(OUTPNG, facecolor="#0f1418", bbox_inches="tight")
print("saved")
