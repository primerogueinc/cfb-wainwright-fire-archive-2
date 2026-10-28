#!/usr/bin/env python3
"""Rubicon Cloud (working model v5.1): annotated YouTube video of the air response to the Oct. 17, 2026 CFB Wainwright fire.

Usage:
  python3 make_video.py <out_dir> [--preview 12.5,40,...]   # --preview writes PNG stills at those video seconds instead
  env: CFB_ARCHIVE_2 (this repository; default: found from this file's location), CFB_ARCHIVE_1 (part 1; default: sibling
       folder ../cfb-wainwright-fire-archive), PR_LOGO (optional Prime Rogue logo PNG)

Inputs (all public, in the two archives):
  repo 2  data/rubicon-ridge/model_v5_1_rubicon_cloud/{context_tracks,wainwright_visits,inbound_response_flights}.csv,
          v5_isochrones_weighted_median.geojson (model), model_v5/front_plane.json (model), model_v5/prefront/*.tif (model),
          satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson (derived vector outline only; no imagery)
  repo 1  data/dem/mrdem30_dtm_3400_10m.tif (NRCan MRDEM-30), data/infrastructure/osm_infrastructure_ODbL.geojson (OSM),
          data/reference/range_standing_orders_2024/rso_points.csv
Output: 1920x1080, 30 fps H.264 (yuv420p, faststart), silent; plus an SRT of the on-screen callouts and a 1280x720 thumbnail.
Frames are rendered with matplotlib + PIL and streamed to ffmpeg over stdin (nothing is written per frame).
"""
import os, sys, json, glob, math, subprocess, argparse, datetime as dt
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon as MPoly
import matplotlib.patheffects as pe
from PIL import Image, ImageDraw, ImageFont
import rasterio, rasterio.features
from rasterio.transform import from_origin
from scipy import ndimage
from pyproj import Transformer

ap = argparse.ArgumentParser()
ap.add_argument("out_dir")
ap.add_argument("--preview", default="")
ap.add_argument("--thumb-only", action="store_true")
ARGS = ap.parse_args()
OUT = ARGS.out_dir
os.makedirs(OUT, exist_ok=True)

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.environ.get("CFB_ARCHIVE_2") or os.path.abspath(os.path.join(HERE, "../../../.."))
if not os.path.isdir(os.path.join(R2, "data/rubicon-ridge")):
    R2 = "/home/claude/cfb-wainwright-fire-archive-2"
R1 = os.environ.get("CFB_ARCHIVE_1") or os.path.join(os.path.dirname(R2), "cfb-wainwright-fire-archive")
LOGO = os.environ.get("PR_LOGO", "")
D51 = f"{R2}/data/rubicon-ridge/model_v5_1_rubicon_cloud"
D5 = f"{R2}/data/rubicon-ridge/model_v5"
sys.path.insert(0, os.path.join(R2, "scripts/rubicon_ridge/model_v5_1"))
from v51_common import PLACES  # noqa: E402  (place coordinates used by v5.1)

W, H, FPS = 1920, 1080, 30
FOCUS = (1300, 450)                       # screen pixel where the camera target sits (free map area right of the legend)
T0 = dt.datetime(2026, 10, 17, 15, 0)     # sim minute 0 = 15:00 MDT Oct. 17

# ---------------------------------------------------------------- palette / fonts
BG, LAND, INK, MUTED = "#0a0f14", "#121a22", "#eef2f5", "#9fb0bf"
PANEL = (20, 28, 37, 255)
ACC = "#ffb347"; FRONT = "#6ec6ff"; BURN_EDGE = "#ff7a4d"
FD = "/usr/share/fonts/opentype/inter"
for f in glob.glob(f"{FD}/Inter-*.otf"):
    font_manager.fontManager.addfont(f)
plt.rcParams["font.family"] = "Inter"
def F(sz, w="Regular", disp=False):
    return ImageFont.truetype(f"{FD}/{'InterDisplay' if disp else 'Inter'}-{w}.otf", sz)
PT = 0.72                                 # px -> pt at dpi 100

TR = Transformer.from_crs(4326, 3400, always_xy=True)
TRI = Transformer.from_crs(3400, 4326, always_xy=True)
def P(lat, lon):
    return TR.transform(lon, lat)
FIRE = (52.7645, -110.9052); FX, FY = P(*FIRE)

# ---------------------------------------------------------------- places
LOCAL = {  # name, lat, lon, marker, label, (dx, dy px), km
    "CYWV": (52.7967, -110.8587, "s", "Airport CYWV\n4.7 km from fire", (16, -64)),
    "CFF7": (52.8280, -110.9045, "H", "Camp helipad (CFF7) · 7 km", (-14, 16)),
    "CWH2": (52.8430, -110.8647, "H", "Health Centre heliport (CWH2) · 9 km", (-14, 16)),
    "CZWN": (52.829, -111.096, "s", "Field 21 (CZWN)\n14.7 km", (-12, -54)),
}
TOWN = (52.8334, -110.8601)
REG = {  # label, key in v51 PLACES, label offset px, ha
    "Edmonton Intl": ("CYEG Edmonton Intl", (-14, -34), "right"),
    "Calgary Intl": ("CYYC Calgary Intl", (14, 14), "left"),
    "Saskatoon": ("CYXE Saskatoon", (-14, -10), "right"),
    "Lloydminster": ("CYLL Lloydminster", (-14, 8), "right"),
    "Lac La Biche": ("CYLB Lac La Biche", (14, 8), "left"),
    "Fort McMurray": ("CYMM Fort McMurray", (-14, 8), "right"),
}
HOSP = {"U of A Hospital site": ("Edmonton, University Hospital", (-16, -2), "right"),
        "Foothills Medical Centre site": ("Calgary, Foothills Medical Centre", (18, -24), "left")}

# ---------------------------------------------------------------- data
fp = json.load(open(f"{D5}/front_plane.json"))
FT0 = fp["fit"]["t_fire_min"]; GB, GC = fp["fit"]["gradient_min_per_km"]
STATIONS = [(s["name"], s["lat"], s["lon"], int(s["t_mdt"][:2]) * 60 + int(s["t_mdt"][3:]) - 900) for s in fp["stations"]]
burn = gpd.read_file(f"{R2}/data/rubicon-ridge/satellogic_2026-10-20/burn_outline_oct17_satellogic.geojson").to_crs(3400)
iso = gpd.read_file(f"{D51}/v5_isochrones_weighted_median.geojson").to_crs(3400)
osm = gpd.read_file(f"{R1}/data/infrastructure/osm_infrastructure_ODbL.geojson").to_crs(3400)
rso = pd.read_csv(f"{R1}/data/reference/range_standing_orders_2024/rso_points.csv")
rso = rso[rso.name.str.match(r"^Range \d+$")]
tracks = pd.read_csv(f"{D51}/context_tracks.csv")
visits = pd.read_csv(f"{D51}/wainwright_visits.csv")

def mins(ts):
    d = pd.to_datetime(ts, format="mixed") - pd.Timestamp(T0)
    return d.dt.total_seconds() / 60.0 if isinstance(d, pd.Series) else d.total_seconds() / 60.0

# ---------------------------------------------------------------- fire arrival field (v5 weighted-median isochrones, model)
def fire_field():
    r0 = rasterio.open(sorted(glob.glob(f"{D5}/prefront/*.tif"))[0])
    res = 20.0
    left, bottom, right, top = r0.bounds
    nx, ny = int((right - left) / res), int((top - bottom) / res)
    tf = from_origin(left, top, res, res)
    inside = rasterio.features.rasterize([burn.geometry.iloc[0]], (ny, nx), transform=tf).astype(bool)
    dom = rasterio.features.rasterize([burn.geometry.iloc[0].buffer(160)], (ny, nx), transform=tf).astype(bool)
    levels = sorted(iso.minutes_after_1500.unique())
    lm = {L: rasterio.features.rasterize(list(iso[iso.minutes_after_1500 == L].geometry), (ny, nx), transform=tf,
                                         all_touched=True).astype(bool) for L in levels}
    anyline = np.zeros((ny, nx), bool)
    for L in levels:
        anyline |= lm[L]
    dist = {L: ndimage.distance_transform_edt(~lm[L]) * res for L in levels}
    lab, n = ndimage.label(dom & ~anyline)
    # pre-front arrival: median of the four v5 pre-front cases (hours MDT) -> minutes after 15:00, resampled to 20 m
    pre = []
    for f in sorted(glob.glob(f"{D5}/prefront/*.tif")):
        a = rasterio.open(f).read(1).astype(float)
        pre.append(np.kron(a, np.ones((2, 2)))[:ny, :nx])
    pre = np.array(pre)
    npre = np.isfinite(pre).sum(axis=0)
    with np.errstate(all="ignore"):
        pre = (np.nanmedian(pre, axis=0) - 15.0) * 60.0
    pre[npre < 2] = np.nan                      # burned before the front in at least two of the four v5 pre-front cases
    T = np.full((ny, nx), np.nan)
    prefront = np.isfinite(pre) & (pre < FT0)
    for k in range(1, n + 1):
        reg = lab == k
        ring = ndimage.binary_dilation(reg, iterations=2) & ~reg
        S = [L for L in levels if (lm[L] & ring).any()]
        if not S:
            continue
        if (prefront & reg).sum() > 0.2 * reg.sum():
            # the region upwind of the first (15:30) contour: the pre-front fire, then front arrival to 15:30
            T[reg] = np.where(prefront[reg], pre[reg], 30.0 - 2.0 * np.clip(dist[levels[0]][reg] / 500.0, 0, 1))
            continue
        if len(S) == 1 and S[0] == levels[0]:
            # upwind of the first (15:30) contour: burned before the front, or between front arrival and 15:30
            t = np.where(np.isfinite(pre[reg]) & (pre[reg] < 29.5), pre[reg], 29.5 - np.minimum(dist[S[0]][reg] / 200.0, 1.5))
            T[reg] = t
        elif len(S) == 1:
            continue                                   # downwind of 16:00: shown by the final outline only
        else:
            ds = np.array([dist[L][reg] for L in S])
            o = np.argsort(ds, axis=0)
            a, b = o[0], o[1]
            La, Lb = np.array(S)[a], np.array(S)[b]
            da = np.take_along_axis(ds, a[None], 0)[0]; db = np.take_along_axis(ds, b[None], 0)[0]
            T[reg] = La + (Lb - La) * da / np.maximum(da + db, 1e-6)
    lines = ndimage.distance_transform_edt(np.isnan(T), return_distances=False, return_indices=True)
    T = np.where(anyline & dom, T[lines[0], lines[1]], T)
    T[~inside] = np.nan
    return T, (left, right, bottom, top)
TF, TF_EXT = fire_field()

# ---------------------------------------------------------------- hillshade
def hillshade():
    with rasterio.open(f"{R1}/data/dem/mrdem30_dtm_3400_10m.tif") as d:
        z = d.read(1, masked=True, out_shape=(d.height // 3, d.width // 3)).astype(float).filled(np.nan)
        b = d.bounds
    z = np.where(np.isfinite(z), z, np.nanmedian(z))
    dy, dx = np.gradient(z, 30.0)
    slope = np.arctan(2.5 * np.hypot(dx, dy)); asp = np.arctan2(-dx, dy)
    az, alt = np.radians(315), np.radians(40)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - asp)
    hs = np.clip(hs, 0, 1)
    rgb = np.zeros(z.shape + (4,))
    base = np.array([0.10, 0.14, 0.17])
    rgb[..., :3] = base[None, None] * 0.6 + hs[..., None] * np.array([0.30, 0.34, 0.36])
    # soft edge so the 29 km tile does not show a hard square
    yy, xx = np.mgrid[0:z.shape[0], 0:z.shape[1]]
    e = np.minimum.reduce([xx, yy, z.shape[1] - 1 - xx, z.shape[0] - 1 - yy]) / (0.12 * z.shape[0])
    rgb[..., 3] = np.clip(e, 0, 1)
    return rgb, (b.left, b.right, b.bottom, b.top)
HS, HS_EXT = hillshade()

# ---------------------------------------------------------------- aircraft
AC = {  # reg: colour, type label, operator label, marker
    "C-FSUG": ("#ffd23f", "King Air 200", "not identified", "o"),
    "C-FLCI": ("#ff5a5f", "H145 heli", "STARS", "D"),
    "C-GSAH": ("#3ddc97", "King Air 200", "Sask. Air Amb.", "o"),
    "C-GLUR": ("#4cc9f0", "King Air 200", "Can West", "o"),
    "C-GWRS": ("#ff9f43", "King Air 200", "not identified", "o"),
    "C-GSAU": ("#b5e655", "King Air 200", "Sask. Air Amb.", "o"),
    "C-GLUD": ("#8aa4ff", "King Air 350", "Can West", "o"),
    "C-FDRP": ("#ff8fc8", "BK117 heli", "STARS", "D"),
    "C-FLNK": ("#c792ff", "H145 heli", "STARS", "D"),
    "C-GLUL": ("#e2b77c", "King Air 200", "Can West", "o"),
    "C-GSAB": ("#9aa5b1", "Citation Mustang", "transit, not a response", "o"),
}
CS = {"C-FSUG": "C-FSUG", "C-FLCI": "STAR03", "C-GSAH": "SLG2", "C-GLUR": "CWA921", "C-GWRS": "C-GWRS", "C-GSAU": "SLG3",
      "C-GLUD": "CWA935", "C-FDRP": "STAR11", "C-FLNK": "STAR01", "C-GLUL": "CWA927", "C-GSAB": "C-GSAB"}
tracks["m"] = mins(tracks.mdt)
tracks["alt"] = pd.to_numeric(tracks.alt_baro, errors="coerce").fillna(0)
tracks["x"], tracks["y"] = TR.transform(tracks.lon.values, tracks.lat.values)
TRK = {r: g.sort_values("m")[["m", "x", "y", "alt"]].to_numpy() for r, g in tracks.groupby("reg")}
FIRST_SEEN = {}

def ac_state(reg, m, trail):
    a = TRK[reg]; t = a[:, 0]
    i = np.searchsorted(t, m)
    pos = None
    if 0 < i < len(t) and t[i] - t[i - 1] <= 5:
        f = (m - t[i - 1]) / max(t[i] - t[i - 1], 1e-9)
        pos = a[i - 1, 1:4] + f * (a[i, 1:4] - a[i - 1, 1:4])
    elif 0 < i and m - t[i - 1] <= 0.25:
        pos = a[i - 1, 1:4]
    j0 = np.searchsorted(t, m - trail)
    seg = a[j0:i]
    if pos is not None:
        seg = np.vstack([seg, [[m, *pos]]]) if len(seg) else np.array([[m, *pos]])
    return pos, seg

# estimated (below-coverage) legs at Wainwright: dotted from last/first heard to the place, around the estimated time
WPL = {"CYWV Wainwright Airport": (52.7967, -110.8587), "CFB Wainwright airfield (CZWN stn)": (52.829, -111.096),
       "Wainwright town": TOWN}
EST = []
for _, v in visits.iterrows():
    reg = v.reg; th = mins(v.last_or_first_heard_mdt); te = mins(v.est_on_ground_mdt)
    a = TRK[reg]; k = np.argmin(abs(a[:, 0] - th)); px, py = a[k, 1], a[k, 2]
    qx, qy = P(*WPL[v.nearest_wainwright_place])
    EST.append(dict(reg=reg, ev=v.event, th=th, te=te, p=(px, py), q=(qx, qy), other=v.other_end,
                    other_t=mins(v.other_end_mdt)))

def dest_city(s):
    if "Edmonton" in s: return "Edmonton"
    if "Calgary" in s: return "Calgary"
    if "Saskatoon" in s: return "Saskatoon"
    return "other"

# ---------------------------------------------------------------- timeline (video seconds <-> sim minutes)
SEGS = []
def seg(kind, dur=None, m0=None, m1=None, rate=None, ch=0, name=""):
    if kind == "run":
        dur = (m1 - m0) / rate
    t0 = SEGS[-1]["t1"] if SEGS else 0.0
    SEGS.append(dict(kind=kind, t0=t0, t1=t0 + dur, m0=m0, m1=m1 if m1 is not None else m0, rate=rate, ch=ch, name=name))
seg("title", 5.0, m0=10, ch=0)
seg("card", 3.0, m0=10, ch=1)
seg("run", m0=10, m1=22, rate=1.0, ch=1)
seg("run", m0=22, m1=26, rate=2.0, ch=1)
seg("run", m0=26, m1=46, rate=2 / 3, ch=1)
seg("run", m0=46, m1=65, rate=1.5, ch=1)
seg("card", 3.0, m0=65, ch=2)
seg("run", m0=65, m1=110, rate=4.0, ch=2)
seg("card", 3.0, m0=110, ch=3)
seg("run", m0=110, m1=114, rate=2.0, ch=3)
seg("run", m0=114, m1=122, rate=2 / 3, ch=3)
seg("run", m0=122, m1=150, rate=2.0, ch=3)
seg("run", m0=150, m1=160, rate=2 / 3, ch=3)
seg("run", m0=160, m1=210, rate=5.0, ch=3)
seg("card", 3.0, m0=210, ch=4)
seg("run", m0=210, m1=780, rate=10.0, ch=4)
seg("hold", 3.0, m0=780, ch=4)
seg("caveats", 8.0, m0=780, ch=5)
seg("credits", 8.0, m0=780, ch=6)
TOTAL = SEGS[-1]["t1"]
NFR = int(round(TOTAL * FPS))
CH_NAMES = {1: "The fire turns", 2: "The quiet", 3: "The response", 4: "The shuttle"}
CH_SUB = {1: "3:10 to 4:05 p.m. · slowed down", 2: "3:34 to 5:20 p.m. · sped up", 3: "4:50 to 6:30 p.m.",
          4: "6:30 p.m. to 4 a.m. · sped up"}

def seg_at(t):
    for s in SEGS:
        if s["t0"] <= t < s["t1"]:
            return s
    return SEGS[-1]

def sim_at(t):
    s = seg_at(t)
    if s["kind"] == "run":
        return s["m0"] + (t - s["t0"]) * s["rate"]
    return s["m0"]

def vt(m, ch=None):
    for s in SEGS:
        if s["kind"] == "run" and s["m0"] <= m < s["m1"] + 1e-9 and (ch is None or s["ch"] == ch):
            return s["t0"] + (m - s["m0"]) / s["rate"]
    raise ValueError(m)

def chap_start(ch):
    return [s for s in SEGS if s["ch"] == ch][0]["t0"]

# ---------------------------------------------------------------- camera (eased; zoom interpolated in log space)
def LL(lat, lon, h):
    x, y = P(lat, lon); return (x, y, h * 1000)
def FO(dxkm, dykm, h):
    return (FX + dxkm * 1000, FY + dykm * 1000, h * 1000)
CAM = [
    (0.0, FO(-5.5, -1.6, 190)),
    (vt(17, 1), FO(-5.5, -1.6, 190)),
    (vt(25.5, 1), FO(-4.9, 4.25, 27)),
    (vt(28.5, 1), FO(-4.9, 4.25, 27)),
    (vt(32, 1), FO(0.1, 0.5, 11)),
    (vt(47, 1), FO(0.1, 0.5, 11)),
    (vt(53, 1), FO(-4.9, 4.25, 27)),
    (chap_start(2) + 3.0, FO(-4.9, 4.25, 27)),
    (vt(78, 2), FO(6, -8, 255)),
    (vt(110, 2) - 0.01, FO(6, -8, 255)),
    (vt(110, 3), FO(6, -8, 255)),
    (vt(116, 3), LL(54.65, -109.6, 760)),
    (vt(126, 3), LL(54.65, -109.6, 760)),
    (vt(144, 3), LL(52.9, -111.05, 170)),
    (vt(159.5, 3), LL(52.9, -111.05, 170)),
    (vt(175, 3), LL(54.65, -109.6, 760)),
    (vt(209.9, 3), LL(54.65, -109.6, 760)),
    (vt(214, 4), LL(52.45, -109.46, 640)),
    (TOTAL, LL(52.45, -109.46, 640)),
]
def ease(u):
    return u * u * (3 - 2 * u) if u < 1 else 1.0
def cam_at(t):
    for (ta, a), (tb, b) in zip(CAM[:-1], CAM[1:]):
        if ta <= t <= tb:
            u = ease((t - ta) / max(tb - ta, 1e-9)) if tb > ta else 1.0
            la, lb = math.log(a[2]), math.log(b[2])
            h = math.exp(la + (lb - la) * u)
            # move the centre in proportion to zoom progress so the motion feels like one continuous flight
            if abs(lb - la) > 1e-6:
                w = (math.log(h) - la) / (lb - la)
            else:
                w = u
            return (a[0] + (b[0] - a[0]) * w, a[1] + (b[1] - a[1]) * w, h)
    return CAM[-1][1]

def limits(c):
    x, y, h = c; s = h / H
    return (x - FOCUS[0] * s, x + (W - FOCUS[0]) * s, y - (H - FOCUS[1]) * s, y + FOCUS[1] * s)

# ---------------------------------------------------------------- callouts (lower thirds + SRT)
def clk(m, ampm=True):
    t = T0 + dt.timedelta(minutes=m)
    return t.strftime("%-I:%M %p").replace("AM", "a.m.").replace("PM", "p.m.")
CALL = [  # (sim minute, chapter, tag, text, seconds on screen)
    (10.2, 1, "3:10 p.m.", "GOES-18 satellite sees the fire: 76 MW of fire power", 5),
    (12.5, 1, "Cold front", "A cold front sweeps in from the northwest at about 77 km/h", 5),
    (16.0, 1, "3:16 p.m.", "Front passes the Field 21 weather station (CZWN): gust 82 km/h", 5),
    (20.0, 1, "3:20 p.m.", "The satellite sees the fire's heat from 3:20 to 3:50 p.m.", 5),
    (27.9, 1, "3:28 p.m.", "Cold front reaches the fire (working model v5)", 5),
    (30.4, 1, "3:30 p.m.", "GOES-18 fire power jumps to 1,215 MW", 5),
    (33.0, 1, "Working model", "Orange area: fire spread, working model v5 estimate (weighted median)", 6),
    (39.0, 1, "3:39 p.m.", "Range 2 / Range 15 (model median 3:39 / 3:42); working-model estimates, 10–90%: 3:34–3:58", 7),
    (41.5, 1, "3:40 p.m.", "GOES-18 fire power: 336 MW", 5),
    (48.0, 1, "Burn outline", "Dashed line: final burn outline, measured from Oct. 20 imagery", 5),
    (56.5, 1, "Oct. 22 briefing", "DND: incident “just before 4:00” (Oct. 22 briefing)", 6),
    (60.5, 1, "Casualties", "Four soldiers died. DND said 39 were injured (as of Oct. 22).", 6),
    (65.3, 2, "3:34–5:20 p.m.", "No aircraft at or below 12,000 ft heard within 80 km of the fire", 5),
    (80.0, 2, "Overhead", "Only airliners passed overhead, at 17,900 ft and higher", 4.5),
    (88.0, 2, "Caveat", "Military aircraft usually don't broadcast; the RCAF helicopters DND mentioned are not visible", 6),
    (114.4, 3, "~4:55 p.m.", "King Air C-FSUG climbs out of the Lac La Biche area toward Wainwright (time estimated)", 5),
    (118.25, 3, "4:58 p.m.", "STARS helicopter C-FLCI (STAR03) lifts off at Edmonton International", 6),
    (124.0, 3, "~40 min", "Flight time Edmonton to Wainwright: about 40 minutes", 5),
    (140.0, 3, "5:20 p.m.", "First aircraft heard within 80 km. The quiet lasted 106 minutes", 5),
    (150.5, 3, "~5:35 p.m.", "C-FSUG and STAR03 reach the Wainwright area (landings estimated, ~5:35–5:37)", 6),
    (154.5, 3, "5:34 p.m.", "SLG2 takes off from Saskatoon toward Wainwright", 5),
    (156.4, 3, "5:36 p.m.", "CWA921 takes off from Edmonton International toward Wainwright", 5),
    (177.5, 3, "~5:57 p.m.", "C-GWRS climbs out of the Fort McMurray area toward Wainwright (estimated)", 5),
    (210.2, 4, "DND, Oct. 22", "“They were cared for first in Wainwright, then flown to hospital across Alberta by STARS, "
                               "Alberta Health Services paramedics, and Royal Canadian Air Force helicopters.”", 7),
    (247.6, 4, "~7:08 p.m.", "First departure from Wainwright: SLG2 toward Saskatoon (lands 7:41 p.m.)", 5),
    (397.5, 4, "~9:38 p.m.", "STAR03 and STAR11 leave; both tracks end at the University of Alberta Hospital site, Edmonton", 6),
    (540.0, 4, "Midnight", "Fixed-wing shuttles continue to Edmonton, Calgary and Saskatoon", 5),
    (749.0, 4, "3:31 a.m.", "STAR01's track ends at the Foothills Medical Centre site, Calgary", 5),
    (772.5, 4, "~3:53 a.m.", "Last departure in the data: SLG3 toward Saskatoon", 4.5),
]
CAPS = []
for m, ch, tag, txt, d in CALL:
    t0 = vt(m, ch)
    if CAPS:
        t0 = max(t0, CAPS[-1]["t0"] + 1.2)
    CAPS.append(dict(t0=t0, t1=t0 + d, tag=tag, txt=txt, ch=ch))
for i, c in enumerate(CAPS):
    if i + 2 < len(CAPS):
        c["t1"] = min(c["t1"], CAPS[i + 2]["t0"])
    nxt = [s["t0"] for s in SEGS if s["kind"] in ("card", "caveats") and s["t0"] > c["t0"]]
    if nxt:
        c["t1"] = min(c["t1"], nxt[0])
CAPS[-1]["t1"] = SEGS[-2]["t0"]          # last caption stays through the hold

# ---------------------------------------------------------------- matplotlib scene
fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=BG)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor(BG); ax.axis("off")
INSET_PX = (1372, 24, 524, 360)          # x, y, w, h (screen px)
axi = fig.add_axes([INSET_PX[0] / W, 1 - (INSET_PX[1] + INSET_PX[3]) / H, INSET_PX[2] / W, INSET_PX[3] / H])
axi.set_facecolor(LAND); axi.set_xticks([]); axi.set_yticks([])
for s in axi.spines.values():
    s.set_color("#5d6d7c"); s.set_linewidth(2)
ZOOMED = []                               # (artist, kind) with alpha control by zoom
def stroke(w=4):
    return [pe.withStroke(linewidth=w, foreground=BG)]

def add_local(a, inset=False):
    hs = a.imshow(HS, extent=HS_EXT, origin="upper", zorder=0, interpolation="bilinear")
    base = osm[osm.name == "CFB Wainwright"]
    arts = []
    for g, kw in [(base, dict(facecolor="#1c2a20", edgecolor="#4f6f55", lw=1.2, alpha=0.55)),
                  (osm[osm.military.fillna("") == "range"], dict(facecolor="#5a4a1a", edgecolor="#a08a3a", lw=0.6, alpha=0.35)),
                  (osm[osm.kind == "water"], dict(facecolor="#1d3a50", edgecolor="none", alpha=0.9))]:
        n0 = len(a.collections); g.plot(ax=a, zorder=1, **kw); arts += a.collections[n0:]
    rd = osm[(osm.kind == "road") & osm.highway.isin(["motorway", "trunk", "primary", "secondary", "tertiary"])]
    n0 = len(a.collections); rd.plot(ax=a, color="#46525e", lw=0.9, zorder=1); arts += a.collections[n0:]
    edge = []
    for part in burn.boundary.explode(index_parts=False):
        xs, ys = part.xy
        ln, = a.plot(xs, ys, color=BURN_EDGE, lw=2.0 if not inset else 1.2, ls=(0, (5, 3)), zorder=4)
        edge.append(ln)
    return hs, arts, edge

hs_main, loc_arts, edge_main = add_local(ax)
for a in loc_arts:
    ZOOMED.append((a, "local", a.get_alpha() or 1.0))
ZOOMED.append((hs_main, "local", 0.9))
for a in edge_main:
    ZOOMED.append((a, "burn", 1.0))
# fire field image (updated per frame)
fire_img = ax.imshow(np.zeros(TF.shape + (4,)), extent=TF_EXT, origin="upper", zorder=3, interpolation="nearest")
iso_lines = []
for _, r in iso.iterrows():
    xs, ys = r.geometry.xy
    ln, = ax.plot(xs, ys, color="#ffe3b0", lw=1.0, alpha=0.0, zorder=4)
    iso_lines.append((r.minutes_after_1500, ln))
# distance rings from the fire centre
th = np.linspace(0, 2 * np.pi, 361)
for rk in (5, 10, 15, 20):
    ln, = ax.plot(FX + rk * 1000 * np.sin(th), FY + rk * 1000 * np.cos(th), color="#ffb38f", lw=1.0, ls=(0, (4, 4)), zorder=2)
    tx = ax.text(FX + rk * 1000 * np.sin(np.radians(200)), FY + rk * 1000 * np.cos(np.radians(200)), f"{rk} km",
                 color="#ffb38f", fontsize=26 * PT, ha="center", va="center", zorder=5, path_effects=stroke())
    ZOOMED += [(ln, "rings", 0.6), (tx, "rings", 1.0)]
ring80, = ax.plot(FX + 80e3 * np.sin(th), FY + 80e3 * np.cos(th), color="#ffb38f", lw=2.0, ls=(0, (6, 5)), zorder=2)
ring80t = ax.text(FX + 80e3 * np.sin(np.radians(225)), FY + 80e3 * np.cos(np.radians(225)), "80 km", color="#ffb38f",
                  fontsize=30 * PT, ha="center", va="center", zorder=5, path_effects=stroke(5), fontweight="semibold")
# airfields / heliports / town / ranges
for k, (la, lo, mk, lab, off) in LOCAL.items():
    x, y = P(la, lo)
    m1, = ax.plot(x, y, mk, color="white", mec="black", ms=11, zorder=7)
    t1 = ax.annotate(lab, (x, y), xytext=off, textcoords="offset pixels", color=INK, fontsize=28 * PT, zorder=7,
                     ha="left" if off[0] > 0 else "right", va="bottom", path_effects=stroke(5), linespacing=1.1)
    ZOOMED += [(m1, "fields", 1.0), (t1, "fields", 1.0)]
hl, hd = 457, np.radians(127); cx, cy = P(52.7967, -110.8587)
rw, = ax.plot([cx - hl * np.sin(hd), cx + hl * np.sin(hd)], [cy - hl * np.cos(hd), cy + hl * np.cos(hd)], color="white",
              lw=5, solid_capstyle="butt", zorder=6)
ZOOMED.append((rw, "fields", 1.0))
tx, ty = P(*TOWN)
tm, = ax.plot(tx, ty, "o", color=INK, ms=9, zorder=7)
tt = ax.annotate("Wainwright", (tx, ty), xytext=(40, -6), textcoords="offset pixels", color=INK, fontsize=32 * PT,
                 fontweight="semibold", zorder=7, path_effects=stroke(5), va="center")
ZOOMED += [(tm, "fields", 1.0), (tt, "fields", 1.0)]
rx, ry = TR.transform(rso.lon.values, rso.lat.values)
rpts, = ax.plot(rx, ry, "D", color="#ffd54a", mec="black", ms=8, zorder=8, ls="")
r2x, r2y = P(52.772, -110.90089)
r2t = ax.annotate("Range 2 / Range 15", (r2x, r2y), xytext=(28, 34), textcoords="offset pixels", color="#ffd54a",
                  fontsize=30 * PT, fontweight="semibold", zorder=9, path_effects=stroke(5),
                  arrowprops=dict(arrowstyle="-", color="#ffd54a", lw=1.5))
ZOOMED += [(rpts, "ranges", 1.0), (r2t, "r2", 1.0)]
firestar, = ax.plot(FX, FY, "*", color="#ff7a4d", mec="white", mew=1.0, ms=20, zorder=9)
firelab = ax.annotate("CFB Wainwright fire", (FX, FY), xytext=(18, 10), textcoords="offset pixels", ha="left", va="bottom",
                      color="#ffb38f", fontsize=30 * PT, fontweight="semibold", zorder=9, path_effects=stroke(5))
ZOOMED += [(firelab, "regional", 1.0)]
# weather stations (front arrival)
st_art = []
for name, la, lo, tm_ in STATIONS:
    if name in ("Wainwright town (PurpleAir)", "CFB Wainwright CZWN"):
        continue
    x, y = P(la, lo)
    short = name.replace(" AGCM", "").replace(" (PurpleAir)", "").replace(" CWCT", "").replace(" CYLL", "")
    mk, = ax.plot(x, y, "o", color=FRONT, mec="white", ms=10, zorder=7)
    off, ha = {"Bellshill": ((-14, 4), "right"), "Castor": ((-14, 4), "right")}.get(short, ((14, 4), "left"))
    lb = ax.annotate(short, (x, y), xytext=off, textcoords="offset pixels", color=INK, fontsize=28 * PT, zorder=7,
                     path_effects=stroke(5), va="bottom", ha=ha)
    st_art.append((tm_, mk, lb, short))
# regional layers
for glon in range(-120, -100):
    gl = np.linspace(48, 61, 60); gx, gy = TR.transform(np.full_like(gl, float(glon)), gl)
    ln, = ax.plot(gx, gy, color="#1a2632", lw=1.0, zorder=0.5); ZOOMED.append((ln, "regional", 1.0))
for glat in range(48, 62):
    go = np.linspace(-121, -99, 80); gx, gy = TR.transform(go, np.full_like(go, float(glat)))
    ln, = ax.plot(gx, gy, color="#1a2632", lw=1.0, zorder=0.5); ZOOMED.append((ln, "regional", 1.0))
lat = np.linspace(48.5, 60.5, 100); bx, by = TR.transform(np.full_like(lat, -110.0), lat)
border, = ax.plot(bx, by, color="#7c8b99", lw=1.6, ls=(0, (8, 4)), zorder=1)
abx, aby = P(52.05, -110.35); skx, sky = P(52.05, -109.65)
abt = ax.text(abx, aby, "ALBERTA", color="#7c8b99", fontsize=28 * PT, ha="right", zorder=2, fontweight="semibold")
skt = ax.text(skx, sky, "SASKATCHEWAN", color="#7c8b99", fontsize=28 * PT, ha="left", zorder=2, fontweight="semibold")
ZOOMED += [(border, "regional", 0.9), (abt, "regional", 1.0), (skt, "regional", 1.0)]
reg_art = []
for lab, (key, off, ha) in REG.items():
    la, lo, _ = PLACES[key]; x, y = P(la, lo)
    mk, = ax.plot(x, y, "s", color=INK, mec="black", ms=10, zorder=6)
    t = ax.annotate(lab, (x, y), xytext=off, textcoords="offset pixels", ha=ha, color=INK, fontsize=30 * PT,
                    fontweight="semibold", zorder=6, path_effects=stroke(5))
    ZOOMED += [(mk, "regional", 1.0), (t, "regional", 1.0)]
    reg_art.append((lab, mk, t))
hosp_art = []
for lab, (key, off, ha) in HOSP.items():
    la, lo, _ = PLACES[key]; x, y = P(la, lo)
    mk, = ax.plot(x, y, "P", color="#ff6b6b", mec="white", ms=13, zorder=6)
    t = ax.annotate(lab, (x, y), xytext=off, textcoords="offset pixels", ha=ha, color="#ffb3b3", fontsize=28 * PT,
                    zorder=6, path_effects=stroke(5))
    hosp_art += [mk, t]
# cold front line and post-front tint
front_ln, = ax.plot([], [], color=FRONT, lw=4, zorder=5, solid_capstyle="round")
front_glow, = ax.plot([], [], color=FRONT, lw=14, alpha=0.18, zorder=5)
front_poly = MPoly(np.zeros((4, 2)), closed=True, facecolor=FRONT, alpha=0.07, edgecolor="none", zorder=2)
ax.add_patch(front_poly)
front_lab = ax.text(0, 0, "", color=FRONT, fontsize=32 * PT, fontweight="bold", zorder=9, path_effects=stroke(6),
                    rotation_mode="anchor", ha="center", va="bottom")
front_arr = ax.annotate("", xy=(0, 0), xytext=(0, 0), zorder=9,
                        arrowprops=dict(arrowstyle="-|>,head_length=0.9,head_width=0.5", lw=5, color=FRONT))
# aircraft artists
ACA = {}
for reg, (col, *_r) in AC.items():
    lc = LineCollection([], linewidths=3.2, zorder=10, capstyle="round"); ax.add_collection(lc)
    mk, = ax.plot([], [], AC[reg][3], color=col, mec="white", mew=1.6, ms=15 if AC[reg][3] == "o" else 13, zorder=12)
    lb = ax.annotate("", (0, 0), xytext=(16, 10), textcoords="offset pixels", color=col, fontsize=28 * PT,
                     fontweight="semibold", zorder=13, path_effects=stroke(5))
    est, = ax.plot([], [], color=col, lw=2.4, ls=(0, (1.5, 3)), zorder=9)
    lci = LineCollection([], linewidths=2.6, zorder=10); axi.add_collection(lci)
    mki, = axi.plot([], [], AC[reg][3], color=col, mec="white", mew=1.2, ms=11, zorder=12)
    esti, = axi.plot([], [], color=col, lw=2.0, ls=(0, (1.5, 3)), zorder=9)
    ACA[reg] = dict(lc=lc, mk=mk, lb=lb, est=est, lci=lci, mki=mki, esti=esti)
pings = ax.scatter([], [], s=[], facecolors="none", edgecolors=[], linewidths=3, zorder=11)
pingsi = axi.scatter([], [], s=[], facecolors="none", edgecolors=[], linewidths=2.5, zorder=11)
# inset static
burn.plot(ax=axi, facecolor="#7a2412", edgecolor=BURN_EDGE, lw=1.2, zorder=3, alpha=0.9)
for rk in (10, 20):
    axi.plot(FX + rk * 1000 * np.sin(th), FY + rk * 1000 * np.cos(th), color="#ffb38f", lw=1.0, ls=(0, (4, 4)), zorder=2, alpha=0.7)
    axi.text(FX + rk * 1000 * np.sin(np.radians(240)), FY + rk * 1000 * np.cos(np.radians(240)), f"{rk} km", color="#ffb38f",
             fontsize=24 * PT, ha="center", va="center", zorder=5, path_effects=stroke(4))
axi.plot([cx - hl * np.sin(hd), cx + hl * np.sin(hd)], [cy - hl * np.cos(hd), cy + hl * np.cos(hd)], color="white", lw=4, zorder=6)
axi.annotate("CYWV", (cx, cy), xytext=(10, -26), textcoords="offset pixels", color=INK, fontsize=26 * PT, zorder=7,
             path_effects=stroke(4), fontweight="semibold")
zx, zy = P(52.829, -111.096)
axi.plot(zx, zy, "s", color="white", mec="black", ms=8, zorder=7)
axi.annotate("CZWN", (zx, zy), xytext=(-10, 10), textcoords="offset pixels", ha="right", color=INK, fontsize=26 * PT,
             zorder=7, path_effects=stroke(4))
axi.plot(tx, ty, "o", color=INK, ms=7, zorder=7)
axi.annotate("Wainwright", (tx, ty), xytext=(10, 8), textcoords="offset pixels", color=INK, fontsize=26 * PT, zorder=7,
             path_effects=stroke(4))
ICX, ICY = P(52.815, -110.90); IW = 52e3; IH = IW * INSET_PX[3] / INSET_PX[2]
axi.set_xlim(ICX - IW / 2, ICX + IW / 2); axi.set_ylim(ICY - IH / 2, ICY + IH / 2)

# ---------------------------------------------------------------- per-frame scene update
def zfade(h, lo, hi):
    """1 when h <= lo, 0 when h >= hi (log-linear)."""
    if h <= lo: return 1.0
    if h >= hi: return 0.0
    return 1 - (math.log(h) - math.log(lo)) / (math.log(hi) - math.log(lo))

def fire_rgba(m):
    out = np.zeros(TF.shape + (4,))
    ok = np.isfinite(TF) & (TF <= m)
    age = np.where(ok, m - TF, 0)
    hot = np.array([1.0, 0.86, 0.35]); mid = np.array([1.0, 0.45, 0.12]); old = np.array([0.55, 0.11, 0.07])
    u = np.clip(age / 4.0, 0, 1)[..., None]; v = np.clip((age - 4) / 10.0, 0, 1)[..., None]
    col = np.where(age[..., None] < 4, hot * (1 - u) + mid * u, mid * (1 - v) + old * v)
    out[..., :3] = col; out[..., 3] = np.where(ok, 0.92, 0)
    return out

def front_geom(m, lim):
    # t = FT0 + GB*dx + GC*dy (min; km from the fire centre)  ->  line where t = m
    g = np.array([GB, GC]); gn = g / np.linalg.norm(g)
    d = (m - FT0) / np.linalg.norm(g)                  # km along gradient (toward later times = downwind)
    p0 = np.array([FX, FY]) + gn * d * 1000
    tang = np.array([-gn[1], gn[0]])
    L = 3000e3
    a, b = p0 - tang * L, p0 + tang * L
    poly = np.array([a, b, b - gn * L, a - gn * L])  # behind the front (already passed)
    return a, b, poly, p0, tang

def update(t, m, ch, cam, kind):
    x0, x1, y0, y1 = limits(cam); h = cam[2] / 1000
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    loc = zfade(h, 40, 110); reg = 1 - zfade(h, 60, 150)
    for a, k, a0 in ZOOMED:
        f = {"local": loc, "burn": max(loc, 0.0) if ch <= 2 else max(loc, 0.6), "rings": zfade(h, 30, 60) if ch in (1, 2) else 0,
             "fields": zfade(h, 30, 45), "ranges": zfade(h, 30, 45), "r2": 1.0 if (ch == 1 and m >= 38.5 and h < 30) else 0.0,
             "regional": reg if ch > 1 or a is border or a is firelab else 0.0}[k]
        a.set_alpha(a0 * f); a.set_visible(f > 0.01)
    firestar.set_visible(True)
    show80 = ch in (2, 3) and m < 150
    f80 = 1 - zfade(h, 120, 230)
    ring80.set_visible(show80 and f80 > 0.01); ring80.set_alpha(0.85 * f80)
    ring80t.set_visible(show80 and f80 > 0.01); ring80t.set_alpha(f80)
    for a in hosp_art:
        a.set_visible(ch >= 4 and reg > 0.2)
    for lab, mk, tt_ in reg_art:
        vis = reg > 0.05 and not (ch == 1) and not (ch == 4 and lab == "Lloydminster")
        px_, py_ = ax.transData.transform(tt_.xy) if vis else (0, 0)
        py_ = H - py_
        ok_ = 50 < py_ < 1040
        mk.set_visible(vis and px_ > 700); tt_.set_visible(vis and ok_ and px_ > (700 if tt_.get_ha() == "left" else 1000))
    if ch >= 3:
        firelab.xyann = (-18, 10); firelab.set_ha("right"); firelab.set_text("Wainwright")
        if ch == 4:
            firelab.set_visible(False)
    else:
        firelab.xyann = (18, 10); firelab.set_ha("left"); firelab.set_text("CFB Wainwright fire")
    # stations: only in chapter 1 at the mid view
    stv = ch == 1 and zfade(h, 120, 160) < 0.99
    sf = 1 - zfade(h, 60, 140)
    for tm_, mk, lb, short in st_art:
        mk.set_visible(stv); lb.set_visible(stv)
        if stv:
            passed = m >= tm_
            mk.set_color(FRONT if passed else "#5d6d7c"); mk.set_alpha(sf); lb.set_alpha(sf)
            lb.set_text(f"{short}\nfront {clk(tm_)}" if passed else short)
    # cold front
    if ch == 1:
        a, b, poly, p0, tang = front_geom(m, (x0, x1, y0, y1))
        front_ln.set_data([a[0], b[0]], [a[1], b[1]]); front_glow.set_data([a[0], b[0]], [a[1], b[1]])
        front_poly.set_xy(poly); front_ln.set_visible(True); front_glow.set_visible(True); front_poly.set_visible(True)
        # label anchored where the line crosses the upper-middle of the free map area
        s = cam[2] / H
        target = np.array([x0 + 1000 * s, y1 - 600 * s])
        q = p0 + tang * np.dot(target - p0, tang)
        ang = math.degrees(math.atan2(tang[1], tang[0]))
        if ang > 90: ang -= 180
        if ang < -90: ang += 180
        inview = x0 + 720 * s < q[0] < x1 - 150 * s and y0 + 250 * s < q[1] < y1 - 80 * s
        front_lab.set_position(q); front_lab.set_rotation(ang); front_lab.set_text("cold front · ~77 km/h")
        front_lab.set_visible(inview and m < 50)
        gn = np.array([GB, GC]) / np.linalg.norm([GB, GC])
        front_arr.xy = tuple(q + gn * 110 * s - tang * 0 * s); front_arr.set_position(tuple(q + gn * 12 * s))
        front_arr.set_visible(inview and m < 50)
    else:
        for a in (front_ln, front_glow, front_poly, front_lab, front_arr):
            a.set_visible(False)
    # fire spread (model)
    if ch <= 2:
        fire_img.set_data(fire_rgba(min(m, 65))); fire_img.set_visible(True); fire_img.set_alpha(max(loc, 0.0))
    else:
        fire_img.set_data(fire_rgba(65)); fire_img.set_visible(loc > 0.01)
    for L, ln in iso_lines:
        ln.set_alpha(0.55 * loc if (ch <= 2 and m >= L) else 0)
    # aircraft
    trail = {1: 8, 2: 15, 3: 30, 4: 50}.get(ch, 30)
    pts, sz, cols = [], [], []
    ptsi, szi, colsi = [], [], []
    active = []
    for reg_, A in ACA.items():
        show = (reg_ == "C-GSAB" and ch == 1) or (reg_ != "C-GSAB" and ch >= 3)
        pos, segp = ac_state(reg_, m, trail) if show else (None, np.zeros((0, 4)))
        col = matplotlib.colors.to_rgb(AC[reg_][0])
        if len(segp) >= 2:
            xy = segp[:, 1:3]; tt_ = segp[:, 0]
            ok = np.diff(tt_) <= 5
            segs = np.stack([xy[:-1], xy[1:]], axis=1)[ok]
            age = (m - tt_[1:])[ok]
            al = np.clip(1 - age / trail, 0, 1) ** 1.3
            A["lc"].set_segments(segs); A["lc"].set_color([(*col, a_) for a_ in al])
            A["lci"].set_segments(segs); A["lci"].set_color([(*col, a_) for a_ in al])
        else:
            A["lc"].set_segments([]); A["lci"].set_segments([])
        if pos is not None:
            A["mk"].set_data([pos[0]], [pos[1]]); A["mk"].set_visible(True)
            A["mki"].set_data([pos[0]], [pos[1]]); A["mki"].set_visible(True)
            alt = "ground" if pos[2] <= 0 else f"{int(round(pos[2] / 25) * 25):,} ft"
            near = math.hypot(pos[0] - cx, pos[1] - cy) < 30e3
            lbl = f"{CS[reg_]} · {alt}"
            if ch == 1:
                lbl = f"{reg_} · {alt}\nbusiness jet in transit"
            A["lb"].xy = (pos[0], pos[1]); A["lb"].set_text(lbl)
            A["lb"].set_visible(not (ch == 4 and near))
            active.append(reg_)
            FIRST_SEEN.setdefault(reg_, t)
        else:
            A["mk"].set_visible(False); A["mki"].set_visible(False); A["lb"].set_visible(False)
        A["est"].set_data([], []); A["esti"].set_data([], [])
    place_labels(ch, kind)
    # estimated legs (below coverage) and pings at Wainwright
    for e in EST if ch >= 3 else []:
        A = ACA[e["reg"]]
        lo_, hi_ = sorted([e["th"], e["te"]])
        if lo_ - 0.01 <= m <= hi_ + 12:
            if e["ev"] == "arrival":
                f = np.clip((m - e["th"]) / max(e["te"] - e["th"], 0.1), 0, 1)
                p, q = np.array(e["p"]), np.array(e["q"])
            else:
                f = np.clip((m - e["te"]) / max(e["th"] - e["te"], 0.1), 0, 1)
                p, q = np.array(e["q"]), np.array(e["p"])
            r = p + (q - p) * f
            A["est"].set_data([p[0], r[0]], [p[1], r[1]]); A["esti"].set_data([p[0], r[0]], [p[1], r[1]])
        col = AC[e["reg"]][0]
        # ping lasts ~1.6 video seconds
        tv = t - vt_safe(e["te"])
        if 0 <= tv <= 1.6:
            pts.append(e["q"]); sz.append((20 + 70 * tv) ** 2); cols.append((*matplotlib.colors.to_rgb(col), 1 - tv / 1.6))
            ptsi.append(e["q"]); szi.append((10 + 40 * tv) ** 2); colsi.append((*matplotlib.colors.to_rgb(col), 1 - tv / 1.6))
        if e["ev"] == "departure":
            tv2 = t - vt_safe(e["other_t"])
            if 0 <= tv2 <= 1.6:
                key = [k for k, v in PLACES.items() if k.split(" ", 1)[-1] in e["other"] or k in e["other"]]
                la_lo = PLACES[key[0]][:2] if key else None
                if la_lo:
                    pts.append(P(*la_lo)); sz.append((20 + 70 * tv2) ** 2)
                    cols.append((*matplotlib.colors.to_rgb(col), 1 - tv2 / 1.6))
    # STAR01 track end at Foothills (not a Wainwright 'departure' with an estimated time)
    tv3 = t - vt_safe(751.2)
    if ch == 4 and 0 <= tv3 <= 1.6:
        pts.append(P(*PLACES["Calgary, Foothills Medical Centre"][:2])); sz.append((20 + 70 * tv3) ** 2)
        cols.append((*matplotlib.colors.to_rgb(AC["C-FLNK"][0]), 1 - tv3 / 1.6))
    if pts:
        pings.set_offsets(np.array(pts)); pings.set_sizes(sz); pings.set_edgecolors(cols); pings.set_visible(True)
    else:
        pings.set_visible(False)
    if ptsi:
        pingsi.set_offsets(np.array(ptsi)); pingsi.set_sizes(szi); pingsi.set_edgecolors(colsi); pingsi.set_visible(True)
    else:
        pingsi.set_visible(False)
    axi.set_visible(ch == 4 and kind in ("run", "hold"))
    return active

def _blocked(ch, kind):
    """Screen boxes (x0, y0, x1, y1; y down) that aircraft labels should avoid: panels and visible place labels."""
    B = [(0, 0, 705, H), (700, 870, W, H)]
    if ch in (2, 3):
        B.append((1465, 0, W, 292))
    if ch == 4:
        B += [(695, 0, 1345, 230), (1365, 0, W, 440)]
    extra = [(None, None, a_) for a_ in hosp_art if hasattr(a_, "xyann")] + [(None, None, lb_) for _, _, lb_, _ in st_art]
    for lab, mk, tt_ in reg_art + [(None, None, firelab)] + extra:
        if tt_.get_visible():
            px_, py_ = ax.transData.transform(tt_.xy); py_ = H - py_
            dx, dy = tt_.xyann; w_ = 0.54 * 30 * max(len(l) for l in tt_.get_text().split("\n"))
            x0_ = px_ + dx if tt_.get_ha() == "left" else px_ + dx - w_
            B.append((x0_ - 2, py_ - dy - 36, x0_ + w_ + 2, py_ - dy + 2))
    return B

def place_labels(ch, kind):
    B = _blocked(ch, kind)
    def hit(b):
        return any(not (b[2] < c[0] or b[0] > c[2] or b[3] < c[1] or b[1] > c[3]) for c in B)
    for reg_, A in ACA.items():
        lb = A["lb"]
        if not lb.get_visible():
            continue
        px_, py_ = ax.transData.transform(lb.xy); py_ = H - py_
        lines = lb.get_text().split("\n")
        w_ = 0.55 * 28 * max(len(l) for l in lines); h_ = 34 * len(lines)
        best = None
        for dx, dy, ha in [(16, 10, "left"), (16, -38 - 36 * (len(lines) - 1), "left"), (-16, 10, "right"),
                           (-16, -38 - 36 * (len(lines) - 1), "right"), (16, 46, "left"), (-16, 46, "right"),
                           (16, -74, "left"), (16, 82, "left")]:
            x0_ = px_ + dx if ha == "left" else px_ + dx - w_
            b = (x0_, py_ - dy - h_, x0_ + w_, py_ - dy)
            if b[0] < 0 or b[2] > W or b[1] < 0:
                continue
            if not hit(b):
                best = (dx, dy, ha, b); break
        if best is None:
            lb.set_visible(False); continue
        lb.xyann = best[:2]; lb.set_ha(best[2]); lb.set_va("bottom")
        B.append(best[3])

_VT = {}
def vt_safe(m):
    if m not in _VT:
        try:
            _VT[m] = vt(m)
        except ValueError:
            _VT[m] = -1e9
    return _VT[m]

# ---------------------------------------------------------------- PIL overlays
def rr(d, box, fill, r=18, outline=None, width=0):
    if fill == PANEL and outline is None:
        outline, width = (52, 66, 82, 255), 2
    d.rounded_rectangle(box, r, fill=fill, outline=outline, width=width)

def wrap(d, txt, font, wmax):
    out, line = [], ""
    for wd in txt.split():
        tst = (line + " " + wd).strip()
        if d.textlength(tst, font=font) <= wmax:
            line = tst
        else:
            out.append(line); line = wd
    out.append(line)
    return out

f_clock = F(92, "Bold", True); f_date = F(32, "Medium"); f_speed = F(29, "Regular"); f_chtag = F(28, "SemiBold")
f_cap = F(38, "Medium"); f_captag = F(30, "Bold"); f_leg = F(28, "Regular"); f_legb = F(28, "Bold"); f_legh = F(28, "SemiBold")
f_cnt_big = F(76, "Bold", True); f_cnt = F(30, "Regular"); f_cnt_h = F(28, "SemiBold")

def speed_text(s):
    if s["kind"] != "run":
        return "Paused"
    r = s["rate"] * 60                       # sim seconds per video second
    if r < 60:
        return f"Speed: 1 second = {r:.0f} seconds"
    v = s["rate"]
    vtxt = f"{v:.0f}" if abs(v - round(v)) < 1e-6 else f"{v:.1f}"
    return f"Speed: 1 second = {vtxt} minute" + ("" if v == 1 else "s")

LOGO_IM = Image.open(LOGO).convert("RGBA") if LOGO and os.path.exists(LOGO) else None

def draw_clock(d, im, t, m, s):
    rr(d, (24, 24, 664, 214), PANEL)
    tt_ = T0 + dt.timedelta(minutes=m)
    hm = tt_.strftime("%-I:%M"); ap_ = tt_.strftime("%p").replace("AM", "a.m.").replace("PM", "p.m.")
    d.text((48, 26), hm, font=f_clock, fill=INK)
    wv = d.textlength(hm, font=f_clock)
    d.text((48 + wv + 14, 60), ap_, font=F(44, "SemiBold", True), fill=INK)
    d.text((48 + wv + 14, 108), "MDT", font=F(28, "Medium"), fill=MUTED)
    d.text((48, 140), tt_.strftime("%a., %b. %-d, %Y").replace("Oct.,", "Oct.").replace("Sat.", "Sat.").replace("Sun.", "Sun."),
           font=f_date, fill=MUTED)
    d.text((48, 178), speed_text(s), font=F(26, "Medium"), fill=ACC)
    if LOGO_IM is not None:
        lg = LOGO_IM.resize((118, 118), Image.LANCZOS); im.alpha_composite(lg, (528, 36))

def draw_chtag(d, ch):
    if ch in CH_NAMES:
        txt = f"{ch}  ·  {CH_NAMES[ch].upper()}"
        w = d.textlength(txt, font=f_chtag)
        rr(d, (24, 228, 24 + w + 40, 280), PANEL, r=12)
        d.text((44, 238), txt, font=f_chtag, fill=ACC)

def draw_legend(d, active_seen, ch, m):
    regs = [r for r in AC if r in active_seen and ch >= 3 and r != "C-GSAB"]
    if not regs:
        return
    y = 296
    hgt = 104 + 44 * len(regs) + (40 if ch >= 3 else 0)
    rr(d, (24, y, 690, y + hgt), PANEL)
    d.text((44, y + 14), "Aircraft (civil transponders)", font=f_legh, fill=INK)
    d.text((44, y + 52), "call sign · reg. · type · operator*", font=F(26, "Regular"), fill=MUTED)
    yy = y + 96
    cw = [max(d.textlength(CS[r] if CS[r] != r else "—", font=f_legb) for r in regs),
          max(d.textlength(r, font=f_leg) for r in regs), max(d.textlength(AC[r][1], font=f_leg) for r in regs)]
    cx_ = [78, 78 + cw[0] + 14, 78 + cw[0] + cw[1] + 28, 78 + cw[0] + cw[1] + cw[2] + 42]
    for r in regs:
        col, typ, op, mk = AC[r]
        if mk == "o":
            d.ellipse((46, yy + 8, 66, yy + 28), fill=col)
        else:
            d.polygon([(56, yy + 5), (68, yy + 18), (56, yy + 31), (44, yy + 18)], fill=col)
        cs = CS[r] if CS[r] != r else "—"
        d.text((cx_[0], yy), cs, font=f_legb, fill=col if cs != "—" else MUTED)
        d.text((cx_[1], yy), r, font=f_leg, fill=INK)
        d.text((cx_[2], yy), typ, font=f_leg, fill=INK)
        d.text((cx_[3], yy), op, font=f_leg, fill=MUTED if op.startswith(("not", "transit")) else INK)
        yy += 44
    if ch >= 3:
        d.text((44, yy + 4), "*inferred from call sign  ·", font=F(26, "Regular"), fill=MUTED)
        xd = 44 + d.textlength("*inferred from call sign  ·", font=F(26, "Regular")) + 14
        for k in range(5):
            d.ellipse((xd + 9 * k, yy + 19, xd + 9 * k + 4, yy + 23), fill=MUTED)
        d.text((xd + 54, yy + 4), "estimated", font=F(26, "Regular"), fill=MUTED)

def draw_caps(d, t):
    vis = [c for c in CAPS if c["t0"] <= t < c["t1"]][-2:]
    if not vis:
        return
    x0, x1 = 700, 1896
    blocks = []
    for c in vis:
        lines = wrap(d, c["txt"], f_cap, x1 - x0 - 60)
        blocks.append((c, lines, 66 + 48 * len(lines)))
    ybot = 1056
    for c, lines, hgt in reversed(blocks):
        a = min(1.0, (t - c["t0"]) / 0.35, (c["t1"] - t) / 0.35)
        y0 = ybot - hgt
        rr(d, (x0, y0, x1, ybot), (10, 15, 20, int(232 * a)))
        d.rectangle((x0, y0 + 14, x0 + 8, ybot - 14), fill=(255, 179, 71, int(255 * a)))
        d.text((x0 + 30, y0 + 14), c["tag"].upper(), font=f_captag, fill=(255, 179, 71, int(255 * a)))
        for i, ln in enumerate(lines):
            d.text((x0 + 30, y0 + 52 + 48 * i), ln, font=f_cap, fill=(238, 242, 245, int(255 * a)))
        ybot = y0 - 12

def draw_goes(d, m):
    x0, y0, x1, y1 = 24, 296, 664, 624
    rr(d, (x0, y0, x1, y1), PANEL)
    d.text((x0 + 22, y0 + 14), "GOES-18 satellite: fire power", font=F(28, "SemiBold"), fill=INK)
    vals = [(10, 76, "3:10"), (20, None, "3:20"), (30, 1215, "3:30"), (40, 336, "3:40"), (50, None, "3:50")]
    bw, gap, base = 90, 30, y1 - 56
    for i, (tm_, v, lab) in enumerate(vals):
        bx = x0 + 34 + i * (bw + gap)
        d.text((bx + bw / 2, base + 10), lab, font=F(26, "Medium"), fill=MUTED, anchor="ma")
        if m >= tm_ + 1:
            if v is None:
                d.text((bx + bw / 2, base - 8), "cloud", font=F(24, "Regular"), fill=MUTED, anchor="md")
            else:
                hh = 160 * v / 1215
                d.rectangle((bx, base - hh, bx + bw, base), fill="#ff7a4d")
                d.text((bx + bw / 2, base - hh - 6), f"{v:,}", font=F(28, "Bold"), fill=INK, anchor="md")
    d.text((x1 - 22, y0 + 52), "MW", font=F(26, "Regular"), fill=MUTED, anchor="ra")

def draw_key(d, m):
    x0, y0 = 24, 640
    items = [("front", "Cold front (v5 fit)"), ("fire", "Fire spread: working model v5 estimate"), ("edge", "Final burn outline (Oct. 20)")]
    rr(d, (x0, y0, 664, y0 + 26 + 46 * len(items)), PANEL)
    yy = y0 + 18
    for k, lab in items:
        if k == "front":
            d.line((x0 + 22, yy + 17, x0 + 66, yy + 17), fill=FRONT, width=6)
        elif k == "fire":
            d.rectangle((x0 + 24, yy + 4, x0 + 64, yy + 30), fill="#ff7a29")
        else:
            for xx in range(x0 + 22, x0 + 66, 14):
                d.line((xx, yy + 17, xx + 8, yy + 17), fill=BURN_EDGE, width=4)
        lines = wrap(d, lab, F(26, "Regular"), 520)
        d.text((x0 + 80, yy + 2), lines[0] if len(lines) == 1 else lab.split(":")[0] + ":", font=F(26, "Regular"), fill=INK)
        if len(lines) > 1:
            d.text((x0 + 80, yy + 32), lab.split(": ")[1], font=F(26, "Regular"), fill=INK)
            yy += 30
        yy += 46

def draw_scale_north(d, cam):
    s = cam[2] / H / 1000                                    # km per px
    target = 420 * s
    nice = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500]
    km = max([n for n in nice if n <= target] or [1])
    px = km / s
    x0 = 60; x1 = x0 + px; y = 1022
    rr(d, (24, 920, 664, 1056), PANEL, r=12)
    d.line((x0, y, x1, y), fill=INK, width=5)
    d.line((x0, y - 12, x0, y + 6), fill=INK, width=4); d.line((x1, y - 12, x1, y + 6), fill=INK, width=4)
    d.text(((x0 + x1) / 2, y - 14), f"{km} km", font=F(28, "SemiBold"), fill=INK, anchor="md")
    # north arrow: true north at the camera centre (grid convergence of EPSG:3400)
    lo, la = TRI.transform(cam[0], cam[1]); xa, ya = P(la + 0.5, lo)
    ang = math.atan2(xa - cam[0], ya - cam[1])
    ncx, ncy = 610, 972
    tip = (ncx + 34 * math.sin(ang), ncy - 34 * math.cos(ang))
    l = (ncx + 16 * math.sin(ang - 2.5), ncy - 16 * math.cos(ang - 2.5))
    r_ = (ncx + 16 * math.sin(ang + 2.5), ncy - 16 * math.cos(ang + 2.5))
    d.polygon([tip, l, (ncx, ncy), r_], fill=INK)
    d.text((ncx - 30, ncy + 4), "N", font=F(28, "Bold"), fill=INK, anchor="rm")

QUIET0, QUIET1 = 34.2, 140.0
def draw_counter(d, m, ch):
    x0, y0, x1 = (1470, 24, 1896) if ch in (2, 3) else (700, 24, 1340)
    if ch == 2 or (ch == 3 and m < 150):
        rr(d, (x0, y0, x1, y0 + 262), PANEL)
        d.text((x0 + 26, y0 + 14), "NO AIRCRAFT HEARD", font=f_cnt_h, fill=ACC)
        d.text((x0 + 26, y0 + 48), "WITHIN 80 KM OF THE FIRE", font=f_cnt_h, fill=ACC)
        big = f"{int(m - QUIET0)} min" if m < QUIET1 else "106 min"
        d.text((x0 + 26, y0 + 82), big, font=f_cnt_big, fill=INK)
        if m < QUIET1:
            d.text((x0 + 26, y0 + 176), "since 3:34 p.m.", font=f_cnt, fill=MUTED)
            d.text((x0 + 26, y0 + 212), "at or below 12,000 ft", font=f_cnt, fill=MUTED)
        else:
            d.text((x0 + 26, y0 + 176), "the quiet ended", font=f_cnt, fill=MUTED)
            d.text((x0 + 26, y0 + 212), "at 5:20 p.m.", font=f_cnt, fill=MUTED)
    elif ch >= 3:
        arr = sum(1 for e in EST if e["ev"] == "arrival" and e["te"] <= m)
        dep = sum(1 for e in EST if e["ev"] == "departure" and e["te"] <= m)
        dest = {"Edmonton": 0, "Calgary": 0, "Saskatoon": 0}
        for e in EST:
            if e["ev"] == "departure" and e["other_t"] <= m:
                dest[dest_city(e["other"])] = dest.get(dest_city(e["other"]), 0) + 1
        fb = F(64, "Bold", True)
        if ch == 3:
            rr(d, (x0, y0, x1, y0 + 262), PANEL)
            d.text((x0 + 26, y0 + 14), "WAINWRIGHT AREA", font=f_cnt_h, fill=ACC)
            d.text((x0 + 26, y0 + 48), "(times estimated)", font=f_cnt_h, fill=ACC)
            for k, (n_, lab_) in enumerate([(arr, "arrivals"), (dep, "departures")]):
                d.text((x0 + 26, y0 + 84 + 86 * k), f"{n_}", font=fb, fill=INK)
                d.text((x0 + 40 + d.textlength(f"{n_}", font=fb), y0 + 112 + 86 * k), lab_, font=F(32, "Medium"), fill=INK)
            return
        fb = F(56, "Bold", True)
        rr(d, (x0, y0, x1, y0 + 200), PANEL)
        d.text((x0 + 26, y0 + 12), "WAINWRIGHT AREA · times estimated", font=f_cnt_h, fill=ACC)
        for k, (n_, lab_) in enumerate([(arr, "arrivals"), (dep, "departures")]):
            d.text((x0 + 26, y0 + 50 + 70 * k), f"{n_}", font=fb, fill=INK)
            d.text((x0 + 106, y0 + 70 + 70 * k), lab_, font=F(32, "Medium"), fill=INK)
        cx2 = x0 + 340
        d.text((cx2, y0 + 54), "Tracks end in:", font=F(28, "Regular"), fill=MUTED)
        for k, c_ in enumerate(["Edmonton", "Calgary", "Saskatoon"]):
            d.text((cx2, y0 + 90 + 34 * k), c_, font=F(28, "Medium"), fill=INK)
            d.text((x1 - 28, y0 + 90 + 34 * k), str(dest[c_]), font=F(28, "Bold"), fill=INK, anchor="ra")

def draw_inset_frame(d):
    x, y, w, h_ = INSET_PX
    rr(d, (x, y + h_ + 4, x + w, y + h_ + 50), PANEL, r=10)
    d.text((x + 16, y + h_ + 12), "Wainwright close-up (52 km wide)", font=F(26, "SemiBold"), fill=INK)

def card_overlay(im, ch, a):
    ov = Image.new("RGBA", (W, H), (6, 9, 12, int(200 * a)))
    d = ImageDraw.Draw(ov)
    col = lambda c, k=1.0: (*c, int(255 * a * k))
    d.text((W / 2, 380), f"CHAPTER {ch}", font=F(40, "SemiBold"), fill=col((255, 179, 71)), anchor="mm")
    d.text((W / 2, 480), CH_NAMES[ch], font=F(110, "Bold", True), fill=col((238, 242, 245)), anchor="mm")
    d.text((W / 2, 590), CH_SUB[ch], font=F(40, "Regular"), fill=col((159, 176, 191)), anchor="mm")
    sub2 = {1: "Satellite fire power, the cold front, and the working model's fire spread",
            2: "What public flight-tracking receivers heard, and didn't",
            3: "The first air ambulances launch and reach Wainwright",
            4: "Aircraft shuttle between Wainwright and hospitals in three cities"}[ch]
    d.text((W / 2, 660), sub2, font=F(36, "Regular"), fill=col((238, 242, 245), 0.9), anchor="mm")
    im.alpha_composite(ov)

def title_frame(t):
    im = Image.new("RGBA", (W, H), BG)
    d = ImageDraw.Draw(im)
    a = min(1, t / 0.6, (5.0 - t) / 0.6)
    if LOGO_IM is not None:
        lg = LOGO_IM.resize((300, 300), Image.LANCZOS)
        lg.putalpha(lg.getchannel("A").point(lambda v: int(v * a)))
        im.alpha_composite(lg, (160, 300))
    c = lambda rgb, k=1.0: (*rgb, int(255 * a * k))
    d.text((520, 300), "Rubicon Cloud", font=F(118, "Bold", True), fill=c((255, 179, 71)))
    d.text((520, 445), "the air response to the", font=F(64, "SemiBold", True), fill=c((238, 242, 245)))
    d.text((520, 525), "CFB Wainwright fire", font=F(64, "SemiBold", True), fill=c((238, 242, 245)))
    d.text((520, 640), "Oct. 17, 2026  ·  public flight data  ·  Prime Rogue Inc.", font=F(38, "Regular"), fill=c((159, 176, 191)))
    d.text((520, 700), "Working model v5.1. All times MDT.", font=F(32, "Regular"), fill=c((159, 176, 191)))
    return im

def text_card(kind, t_in, dur):
    im = Image.new("RGBA", (W, H), BG); d = ImageDraw.Draw(im)
    a = min(1, t_in / 0.5, (dur - t_in) / 0.5)
    c = lambda rgb, k=1.0: (*rgb, int(255 * a * k))
    if kind == "caveats":
        d.text((160, 120), "What this data can't show", font=F(76, "Bold", True), fill=c((255, 179, 71)))
        items = ["Receivers hear nothing below about 4,250 ft (about 2,000 ft above the ground), so landings and "
                 "low hops are invisible.",
                 "Military aircraft don't broadcast, so the RCAF helicopters DND mentioned are not visible.",
                 "Landing and take-off times near Wainwright are estimates (±3–5 min).",
                 "Operators are inferred from public call-sign prefixes. Nothing here says who was aboard any aircraft.",
                 "Fire timing is working model v5, an estimate, not a final finding."]
        y = 270
        for it in items:
            lines = wrap(d, it, F(42, "Regular"), 1500)
            d.ellipse((168, y + 18, 186, y + 36), fill=c((255, 179, 71)))
            for i, ln in enumerate(lines):
                d.text((220, y + 54 * i), ln, font=F(42, "Regular"), fill=c((238, 242, 245)))
            y += 54 * len(lines) + 34
    else:
        d.text((160, 90), "Sources and credits", font=F(70, "Bold", True), fill=c((255, 179, 71)))
        items = ["Aircraft positions: adsb.lol open archive, © adsb.lol contributors (ODbL)",
                 "Fire power: NOAA GOES-18 · Weather: Environment and Climate Change Canada",
                 "Terrain: NRCan MRDEM-30 (Open Government Licence – Canada)",
                 "Map data: © OpenStreetMap contributors (ODbL) · Range grids: 2024 Range Standing Orders",
                 "Statements: DND, Oct. 22, 2026 briefing",
                 "Analysis: Prime Rogue working model v5 / v5.1 (CC BY-NC 4.0)"]
        y = 210
        for it in items:
            d.text((160, y), it, font=F(36, "Regular"), fill=c((238, 242, 245)))
            y += 64
        d.text((160, y + 30), "Data, code and method:", font=F(36, "SemiBold"), fill=c((159, 176, 191)))
        d.text((160, y + 80), "github.com/primerogueinc/cfb-wainwright-fire-archive-2", font=F(44, "Bold"), fill=c((110, 198, 255)))
        if LOGO_IM is not None:
            lg = LOGO_IM.resize((170, 170), Image.LANCZOS)
            lg.putalpha(lg.getchannel("A").point(lambda v: int(v * a)))
            im.alpha_composite(lg, (1690, 860))
        d.text((160, 990), "Prime Rogue Inc. · Calgary", font=F(34, "SemiBold"), fill=c((238, 242, 245)))
    return im

# ---------------------------------------------------------------- frame
SEEN = set()
def render(t):
    s = seg_at(t)
    if s["kind"] == "title":
        return title_frame(t - s["t0"]).convert("RGB")
    if s["kind"] in ("caveats", "credits"):
        return text_card(s["kind"], t - s["t0"], s["t1"] - s["t0"]).convert("RGB")
    m = sim_at(t); ch = s["ch"]; cam = cam_at(t)
    active = update(t, m, ch, cam, s["kind"])
    SEEN.update(active)
    fig.canvas.draw()
    im = Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba(), "raw", "RGBA", 0, 1).copy()
    ov = Image.new("RGBA", (W, H), (0, 0, 0, 0)); d = ImageDraw.Draw(ov)
    draw_clock(d, ov, t, m, s)
    draw_chtag(d, ch)
    seen = {r for r in SEEN if r in AC}
    draw_legend(d, seen, ch, m)
    if ch == 1 and s["kind"] == "run":
        draw_goes(d, m)
        if m >= 27:
            draw_key(d, m)
    if s["kind"] in ("run", "hold"):
        draw_counter(d, m, ch) if (ch >= 2) else None
        if ch == 4:
            draw_inset_frame(d)
        draw_scale_north(d, cam)
        draw_caps(d, t)
    im.alpha_composite(ov)
    if s["kind"] == "card":
        u = t - s["t0"]; a = min(1, u / 0.4, (3.0 - u) / 0.4)
        card_overlay(im, ch, max(a, 0))
    if t < 5.6:                                           # cross-fade from title
        pass
    return im.convert("RGB")

# ---------------------------------------------------------------- SRT
def srt_time(x):
    ms = int(round(x * 1000)); h_, ms = divmod(ms, 3600000); m_, ms = divmod(ms, 60000); s_, ms = divmod(ms, 1000)
    return f"{h_:02d}:{m_:02d}:{s_:02d},{ms:03d}"
def write_srt(path):
    ev = [(0.3, 4.8, "Rubicon Cloud: the air response to the CFB Wainwright fire\nOct. 17, 2026 · public flight data · Prime Rogue Inc.")]
    for s in SEGS:
        if s["kind"] == "card":
            ev.append((s["t0"] + 0.2, s["t1"] - 0.1, f"Chapter {s['ch']}: {CH_NAMES[s['ch']]}\n{CH_SUB[s['ch']]}"))
        if s["kind"] == "caveats":
            ev.append((s["t0"] + 0.2, s["t1"] - 0.1, "What this data can't show: receivers hear nothing below about 4,250 ft; "
                       "military aircraft don't broadcast; landing times are estimates; operators are inferred from call signs."))
        if s["kind"] == "credits":
            ev.append((s["t0"] + 0.2, s["t1"] - 0.1, "Sources: adsb.lol (ODbL), NOAA GOES, ECCC, NRCan MRDEM-30, © OpenStreetMap "
                       "contributors, 2024 Range Standing Orders, Prime Rogue model v5/v5.1 (CC BY-NC 4.0).\n"
                       "github.com/primerogueinc/cfb-wainwright-fire-archive-2"))
    for c in CAPS:
        timed = c["tag"][0].isdigit() or c["tag"].startswith("~")
        ev.append((c["t0"], c["t1"], f"{c['tag']}: {c['txt']}" if timed else c["txt"]))
    ev.sort()
    with open(path, "w") as f:
        for i, (a, b, txt) in enumerate(ev, 1):
            f.write(f"{i}\n{srt_time(a)} --> {srt_time(b)}\n{txt}\n\n")

# ---------------------------------------------------------------- thumbnail
def thumbnail(path):
    # regional map with every response track, rendered from the same scene at the chapter-4 camera, then scaled
    t = vt(779.9, 4); m = 779.9
    cam = (*P(52.9, -112.3), 700e3)
    update(t, m, 4, cam, "thumb")
    for reg_, A in ACA.items():
        if reg_ == "C-GSAB":
            continue
        a = TRK[reg_]; k = (a[:, 0] >= 115) & (a[:, 0] <= 780)
        xy = a[k][:, 1:3]; tt_ = a[k][:, 0]
        ok = np.diff(tt_) <= 5
        A["lc"].set_segments(np.stack([xy[:-1], xy[1:]], axis=1)[ok]); A["lc"].set_color(AC[reg_][0]); A["lc"].set_alpha(1.0)
        A["lc"].set_linewidth(5)
        A["mk"].set_visible(False); A["lb"].set_visible(False); A["est"].set_data([], [])
    pings.set_visible(False); axi.set_visible(False)
    for _, mk, t_ in reg_art:
        t_.set_visible(False)
    for a_ in hosp_art:
        a_.set_visible(False)
    abt.set_visible(False); skt.set_visible(False)
    fig.canvas.draw()
    im = Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba(), "raw", "RGBA", 0, 1).copy()
    d = ImageDraw.Draw(im)
    grad = Image.linear_gradient("L").rotate(-90).resize((W, H))
    shade = Image.new("RGBA", (W, H), (6, 9, 12, 255)); shade.putalpha(grad.point(lambda v: int(min(235, max(0, v * 1.6 - 60)))))
    im.alpha_composite(shade)
    d = ImageDraw.Draw(im)
    d.text((70, 120), "THE AIR", font=F(210, "Black", True), fill=INK)
    d.text((70, 330), "RESPONSE", font=F(210, "Black", True), fill=ACC)
    d.text((76, 590), "CFB Wainwright fire · Oct. 17, 2026", font=F(64, "Bold", True), fill=INK)
    d.text((76, 680), "106 minutes of empty sky, then a night-long shuttle", font=F(50, "SemiBold"), fill=MUTED)
    d.text((76, 880), "Public flight data · working model v5.1", font=F(46, "Medium"), fill=MUTED)
    if LOGO_IM is not None:
        im.alpha_composite(LOGO_IM.resize((250, 250), Image.LANCZOS), (W - 300, H - 300))
    im.convert("RGB").resize((1280, 720), Image.LANCZOS).save(path, optimize=True)

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    print(f"total {TOTAL:.2f} s, {NFR} frames")
    for s in SEGS:
        print(f"  {s['kind']:8s} ch{s['ch']} {s['t0']:7.2f}-{s['t1']:7.2f}  sim {s['m0']}-{s['m1']}")
    write_srt(os.path.join(OUT, "rubicon_cloud_air_response.srt"))
    json.dump(dict(total_s=TOTAL, chapters={str(s["ch"]): round([x for x in SEGS if x["ch"] == s["ch"]][0]["t0"], 2)
                                            for s in SEGS}),
              open(os.path.join(OUT, "timeline.json"), "w"), indent=1)
    if ARGS.preview:
        for tt in [float(x) for x in ARGS.preview.split(",")]:
            # warm SEEN so legends reflect what was visible up to that point
            for t2 in np.arange(0, tt, 0.5):
                s2 = seg_at(t2)
                if s2["kind"] == "run":
                    for reg_ in AC:
                        show = (reg_ == "C-GSAB" and s2["ch"] == 1) or (reg_ != "C-GSAB" and s2["ch"] >= 3)
                        if show and ac_state(reg_, sim_at(t2), 1)[0] is not None:
                            SEEN.add(reg_)
            render(tt).save(os.path.join(OUT, f"preview_{tt:06.1f}.png"))
        sys.exit(0)
    if ARGS.thumb_only:
        thumbnail(os.path.join(OUT, "rubicon_cloud_thumbnail.png")); sys.exit(0)
    mp4 = os.path.join(OUT, "rubicon_cloud_air_response.mp4")
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS),
           "-i", "-", "-f", "lavfi", "-i", "anullsrc=channel_layout=stereo:sample_rate=48000",
           "-c:v", "libx264", "-preset", "slow", "-crf", "22", "-tune", "animation", "-pix_fmt", "yuv420p",
           "-c:a", "aac", "-b:a", "32k", "-shortest", "-movflags", "+faststart", mp4]
    pr = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    import time
    t_start = time.time()
    for i in range(NFR):
        pr.stdin.write(render(i / FPS).tobytes())
        if i % 300 == 0:
            print(f"frame {i}/{NFR}  {time.time() - t_start:.0f}s", flush=True)
    pr.stdin.close(); pr.wait()
    thumbnail(os.path.join(OUT, "rubicon_cloud_thumbnail.png"))
    print("done", mp4)
