"""Rubicon Ridge (working model v5): shared inputs, grid, FBP grass functions and the spread engine.

Reuses the v4 post-front set-up: the model grid, water mask, starting outline (12:17 Landsat + 2:49 p.m. VIIRS hull) and GOES-18
pixel footprints all come from scripts/authentic_stoic_alternative_hypothesis/spread.py (executed up to its "runs" section,
exactly as model_v4/post_front.py and fit_post.py do), so the v5 grid is identical to v4's (40 m, EPSG:3400, meta.npy).

What is new here:
  * spread_v5(): the v4 Huygens/Dijkstra spread (16-neighbour stencil, ellipse with the ignition point at the rear focus),
    compiled with numba, with (a) a front that reaches each cell at its own time tf[i, j] (moving line), (b) a post-front wind
    and head-spread schedule indexed by time since the front reached the fire, and (c) a per-cell spread factor (spring burn).
  * FBP grass: O-1a (matted) and O-1b (standing), ST-X-3 eq. 35a/35b with the 2009 curing function; ISI and hourly FFMC.
  * load_scar(path): the burn scar used for fitting is a parameter (default: data/infrastructure/scar_patches_revised.geojson,
    largest patch, as v4), so a sharper outline (e.g. the Satellogic 70 cm outline) can be swapped in without touching the grid.
"""
import os, math, json
import numpy as np
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.environ.get("V5_SCRATCH", "/tmp/claude-0/-home-claude/11d4f8cb-573a-58d8-9723-d23e8fe40496/scratchpad/v5"), "numba_cache"))
import numba as nb
import heapq

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "../../.."))
# v1-v4 inputs live in part 1 of the archive (primerogueinc/cfb-wainwright-fire-archive).
# Set CFB_ARCHIVE_1 to its checkout; default assumes it sits next to this repository.
V4REPO = os.path.abspath(os.environ.get("CFB_ARCHIVE_1", os.path.join(REPO, "../cfb-wainwright-fire-archive")))
OUT = os.path.join(REPO, "data/rubicon-ridge/model_v5")
SCRATCH = os.environ.get("V5_SCRATCH", "/tmp/claude-0/-home-claude/11d4f8cb-573a-58d8-9723-d23e8fe40496/scratchpad/v5")
DEFAULT_SCAR = os.path.join(V4REPO, "data/infrastructure/scar_patches_revised.geojson")
SPRING_BURN = os.path.join(V4REPO, "data/fuel/burn_2026-04-20_pda2.geojson")
FIRE_LONLAT = (-110.9052, 52.7645)     # fire centre used by the Rubicon Ridge data pull (front-plane origin)

# Reference points (lon, lat). Range grids: 2024 Range Standing Orders (rso_points.csv); road head as given for Rubicon Ridge.
def ref_points():
    import pandas as pd
    r = pd.read_csv(os.path.join(V4REPO, "data/reference/range_standing_orders_2024/rso_points.csv"))
    pick = lambda n: tuple(r.loc[r.name == n, ["lon", "lat"]].iloc[0])
    return {"Range 2 (RSO)": pick("Range 2"), "Range 11 (RSO)": pick("Range 11"), "Range 15 (RSO)": pick("Range 15"),
            "Road head": (-110.89965, 52.77045),
            "Range 2 point (OSM, v4)": (-110.9063, 52.7728)}   # v4's point, kept only for the v4 comparison

# ---------------------------------------------------------------- grid and v4 inputs (exec of spread.py, as in v4)
def load_v4_setup():
    src = open(os.path.join(V4REPO, "scripts/authentic_stoic_alternative_hypothesis/spread.py")).read().split("# ---------------- runs")[0]
    src = src.replace('OUT = os.path.join(REPO, "data/authentic-stoic-alternative-hypothesis")', f'OUT = "{SCRATCH}"')
    g = {"__name__": "spread_v4_setup"}
    import sys
    argv = sys.argv; sys.argv = ["x", V4REPO]
    try:
        exec(src, g)
    finally:
        sys.argv = argv
    keep = ["minx", "maxy", "RES", "H", "W", "T", "CRS", "scar", "hull", "scar_m", "hull_m", "water", "foot", "rast", "mdt"]
    return {k: g[k] for k in keep}

def load_scar(path, S, mode="largest"):
    """Rasterise the fitting scar onto the fixed v4 grid. mode 'largest' (v4 behaviour) or 'union' (all polygons)."""
    import geopandas as gpd
    from shapely.ops import unary_union
    gdf = gpd.read_file(path).to_crs(S["CRS"])
    geom = max(gdf.geometry, key=lambda g: g.area) if mode == "largest" else unary_union(list(gdf.geometry))
    return geom, S["rast"](geom)

def cell_of(S, lon, lat):
    from pyproj import Transformer
    x, y = Transformer.from_crs(4326, S["CRS"], always_xy=True).transform(lon, lat)
    return int((S["maxy"] - y) // S["RES"]), int((x - S["minx"]) // S["RES"]), x, y

def grid_xy(S):
    yy, xx = np.mgrid[0:S["H"], 0:S["W"]]
    return S["minx"] + (xx + 0.5) * S["RES"], S["maxy"] - (yy + 0.5) * S["RES"]

# ---------------------------------------------------------------- FBP grass (ST-X-3; Wotton et al. 2009 curing function)
def curing_factor(c):
    return 0.005 * (math.exp(0.061 * c) - 1) if c < 58.8 else 0.176 + 0.02 * (c - 58.8)

def isi(ffmc, ws):
    fm = 147.27723 * (101 - ffmc) / (59.5 + ffmc)
    sf = 19.115 * math.exp(-0.1386 * fm) * (1 + fm ** 5.31 / 4.93e7)
    fw = math.exp(0.05039 * ws) if ws < 40 else 12 * (1 - math.exp(-0.0818 * (ws - 28)))
    return fw * sf          # 19.115 = 0.208 x 91.9 already: no second 0.208 (the bug fixed in v4)

def ros_grass(I, fuel, curing=95.0):
    """Head rate of spread, m/min. O-1a matted grass: a=190, b=0.0310, c=1.4. O-1b standing grass: a=250, b=0.0350, c=1.7."""
    a, b, c = {"O-1a": (190, 0.0310, 1.4), "O-1b": (250, 0.0350, 1.7)}[fuel]
    return a * (1 - math.exp(-b * I)) ** c * curing_factor(curing)

def hffmc(f0, temp, rh, ws, hours=1.0):     # Van Wagner (1977) hourly FFMC, as part8/ignition_backcalc.py
    mo = 147.27723 * (101 - f0) / (59.5 + f0)
    ed = 0.942 * rh ** 0.679 + 11 * math.exp((rh - 100) / 10) + 0.18 * (21.1 - temp) * (1 - math.exp(-0.115 * rh))
    ew = 0.618 * rh ** 0.753 + 10 * math.exp((rh - 100) / 10) + 0.18 * (21.1 - temp) * (1 - math.exp(-0.115 * rh))
    if mo > ed:
        a = rh / 100; k0 = 0.424 * (1 - a ** 1.7) + 0.0694 * ws ** 0.5 * (1 - a ** 8); k = k0 * 0.0579 * math.exp(0.0365 * temp)
        m = ed + (mo - ed) * 10 ** (-k * hours)
    elif mo < ew:
        a = (100 - rh) / 100; k0 = 0.424 * (1 - a ** 1.7) + 0.0694 * ws ** 0.5 * (1 - a ** 8); k = k0 * 0.0579 * math.exp(0.0365 * temp)
        m = ew - (ew - mo) * 10 ** (-k * hours)
    else:
        m = mo
    return 59.5 * (250 - m) / (147.27723 + m)

def lb_grass(ws):
    return 1.1 * max(ws, 1.0) ** 0.464     # FBP grass length-to-breadth (ST-X-3 eq. 80)

# ---------------------------------------------------------------- spread engine
OFF = [(dy, dx) for dy in range(-2, 3) for dx in range(-2, 3) if (dy, dx) != (0, 0) and math.gcd(abs(dy), abs(dx)) == 1]
def offsets(res):
    dy = np.array([o[0] for o in OFF], np.int64); dx = np.array([o[1] for o in OFF], np.int64)
    dist = np.array([math.hypot(a, b) * res for a, b in OFF]); brg = np.array([math.atan2(b, -a) for a, b in OFF])
    return dy, dx, dist, brg

STEP = 0.25   # schedule sampling, minutes

@nb.njit(cache=True)
def _hpush(hk, hv, n, key, val):
    i = n; hk[i] = key; hv[i] = val
    while i > 0:
        p = (i - 1) // 2
        if hk[p] <= hk[i]:
            break
        hk[p], hk[i] = hk[i], hk[p]; hv[p], hv[i] = hv[i], hv[p]; i = p
    return n + 1

@nb.njit(cache=True)
def _hpop(hk, hv, n):
    n -= 1; hk[0] = hk[n]; hv[0] = hv[n]; i = 0
    while True:
        l = 2 * i + 1; r = l + 1; m = i
        if l < n and hk[l] < hk[m]: m = l
        if r < n and hk[r] < hk[m]: m = r
        if m == i:
            break
        hk[m], hk[i] = hk[i], hk[m]; hv[m], hv[i] = hv[i], hv[m]; i = m
    return n

@nb.njit(cache=True)
def spread_v5(seed_t, tf, fuel, rfac, head_rad, ros, ecc, t_fire, t_end, dy, dx, dist, brg):
    """Time-dependent shortest arrival (minutes after 15:00 MDT).
    seed_t : arrival of the pre-front fire (inf elsewhere). A burning cell spreads post-front only once the front has reached
             it (departure time max(t, tf[i, j])); before that its growth is the pre-front model's (already in seed_t).
    head_rad, ros, ecc : schedule sampled every STEP minutes from the moment the front reaches the fire (t_fire); a cell uses the
             value at (departure - t_fire), so the wind pattern travels with the front.
    rfac   : per-cell spread factor of the cell being entered (spring-burn regrowth)."""
    H, W = seed_t.shape
    # float64 always. The pre-front rasters are float32; a float32 arrival array rounds each float64 candidate time on storage,
    # so "nt < arr" can succeed repeatedly for the same stored value: endless re-pushes (gigabytes of heap) and wrong times.
    # This was the cause of the first v5 attempt's out-of-memory failure and of its (discarded) first block of results.
    arr = seed_t.astype(np.float64)
    # Binary min-heap on preallocated arrays (keys = time, vals = flat cell index); replaces the numba heapq list.
    cap = H * W * 4
    hk = np.empty(cap, np.float64); hv = np.empty(cap, np.int64); n_h = 0
    for i in range(H):
        for j in range(W):
            if np.isfinite(arr[i, j]):
                n_h = _hpush(hk, hv, n_h, arr[i, j], i * W + j)
    n = ros.shape[0]; K = dy.shape[0]
    while n_h > 0:
        t = hk[0]; c = hv[0]
        n_h = _hpop(hk, hv, n_h)
        i = c // W; j = c % W
        if t > arr[i, j]:
            continue
        td = t if t > tf[i, j] else tf[i, j]
        if td > t_end:
            continue
        k = int((td - t_fire) / 0.25)
        if k < 0: k = 0
        if k >= n: k = n - 1
        h = head_rad[k]; R = ros[k]; e = ecc[k]
        for q in range(K):
            a = i + dy[q]; b = j + dx[q]
            if a < 0 or a >= H or b < 0 or b >= W or not fuel[a, b]:
                continue
            r = R * rfac[a, b] * (1 - e) / (1 - e * math.cos(brg[q] - h))
            nt = td + dist[q] / r
            if nt < arr[a, b]:
                arr[a, b] = nt
                if n_h >= hk.shape[0]:
                    hk2 = np.empty(2 * hk.shape[0], np.float64); hv2 = np.empty(2 * hk.shape[0], np.int64)
                    hk2[:n_h] = hk[:n_h]; hv2[:n_h] = hv[:n_h]; hk = hk2; hv = hv2
                n_h = _hpush(hk, hv, n_h, nt, a * W + b)
    return arr

def hm(t):
    if t is None or not np.isfinite(t): return None
    m = int(round(15 * 60 + float(t))); return f"{m // 60}:{m % 60:02d}"

def wquant(x, w, q):
    x = np.asarray(x, float); w = np.asarray(w, float); ok = np.isfinite(x)
    if not ok.any(): return None
    x, w = x[ok], w[ok]; o = np.argsort(x); x, w = x[o], w[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q, c, x))
