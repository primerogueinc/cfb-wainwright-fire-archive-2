"""Rubicon Cloud (working model v5.1): shared paths, places and ADS-B helpers.

v5.1 sits alongside Rubicon Ridge (model v5) and does not change any v5 number. It adds aircraft positions (ADS-B) to the
timeline. Inputs:
  * data/rubicon-ridge/raw/adsb_night_low/wainwright_fire_night_low_near_base_all_positions.csv  (supplied by Prime Rogue,
    1,671 positions; shown by v51_air.py to be a filtered subset of the adsb.lol archive already kept in part 1)
  * part 1: data/adsb/oct17/raw_traces/  (adsb.lol open archive, ODbL 1.0, original readsb trace files) for the full flights
  * model v5: data/rubicon-ridge/model_v5/ (front plane, run weights, arrival times) and, for the isochrone map only, the stored
    scar-limited arrival grids (V5_SCRATCH/grids; not in the repository)
Times are MDT (UTC-6) unless a column says utc.
"""
import os, sys, json, glob, gzip, math
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "model_v5"))
import v5_common as C5                                   # noqa: E402  (REPO, V4REPO, FIRE_LONLAT, hm, wquant, SCRATCH)

REPO, V4REPO = C5.REPO, C5.V4REPO
OUT = os.path.join(REPO, "data/rubicon-ridge/model_v5_1_rubicon_cloud")
V5 = os.path.join(REPO, "data/rubicon-ridge/model_v5")
SUPPLIED = os.path.join(REPO, "data/rubicon-ridge/raw/adsb_night_low/wainwright_fire_night_low_near_base_all_positions.csv")
TRACES = os.path.join(V4REPO, "data/adsb/oct17/raw_traces")
BURN = os.path.join(REPO, "data/rubicon-ridge/satellogic_2026-10-20/outline_for_model_v5.geojson")   # vectors only
FIRE = (52.7645, -110.9052)              # fire centre (lat, lon): the v5 front-plane origin
FRONT_AT_FIRE = 27.89                    # minutes after 15:00 MDT (v5 front_plane.json; leave-one-out 27.5-28.5)
MDT_OFF = pd.Timedelta(hours=6)

# Places used to label where tracks start and end. Coordinates are approximate (published aerodrome reference points rounded to
# ~0.001 deg; hospital helipads to the hospital site). They are used only for "nearest place" labels and distances.
PLACES = {
    "CYWV Wainwright Airport":                (52.7967, -110.8587, 2230),
    "CFB Wainwright airfield (CZWN stn)":     (52.8294, -111.0956, 2260),
    "Wainwright town":                        (52.8390, -110.8573, 2200),
    "CYEG Edmonton Intl":                     (53.3097, -113.5797, 2373),
    "CYYC Calgary Intl":                      (51.1225, -114.0133, 3606),
    "CYXE Saskatoon":                         (52.1708, -106.6997, 1653),
    "CYLL Lloydminster":                      (53.3092, -110.0725, 2193),
    "CYMM Fort McMurray":                     (56.6533, -111.2219, 1211),
    "CYXH Medicine Hat":                      (50.0189, -110.7208, 2352),
    "CYOD Cold Lake":                         (54.4050, -110.2794, 1775),
    "CYQF Red Deer":                          (52.1822, -113.8944, 2968),
    "CYQR Regina":                            (50.4319, -104.6658, 1894),
    "CYPA Prince Albert":                     (53.2142, -105.6728, 1405),
    "CYLB Lac La Biche":                      (54.7703, -112.0317, 1888),
    "CYBF Bonnyville":                        (54.3042, -110.7444, 1838),
    "CZSP St. Paul":                          (53.9933, -111.3800, 2147),
    "Edmonton, University Hospital":          (53.5206, -113.5236, 2200),
    "Edmonton, Royal Alexandra Hospital":     (53.5577, -113.4970, 2200),
    "Calgary, Foothills Medical Centre":      (51.0650, -114.1340, 3500),
}
WAINWRIGHT_KEYS = ["CYWV Wainwright Airport", "CFB Wainwright airfield (CZWN stn)", "Wainwright town"]

# Operators from public ICAO designators / call-sign prefixes (Nav Canada list of Canadian ICAO designators; Wikipedia list of
# Canadian airlines). Mapping a call sign to an operator is an INFERENCE from the prefix; nothing here says who was aboard.
CALLSIGN_OPERATORS = {
    "STAR": ("STARS (Shock Trauma Air Rescue Service), air ambulance helicopters", "STAR"),
    "CWA": ("Can West Corporate Air Charters (call sign CAN WEST); flies fixed-wing air ambulance under contract to "
            "Alberta Health Services", "CWA"),
    "SLG": ("Saskatchewan Government Air Ambulance Service (call sign LIFEGUARD)", "SLG"),
}


def hav(la1, lo1, la2, lo2):
    p = np.pi / 180
    a = 0.5 - np.cos((la2 - la1) * p) / 2 + np.cos(la1 * p) * np.cos(la2 * p) * (1 - np.cos((lo2 - lo1) * p)) / 2
    return 12742 * np.arcsin(np.sqrt(a))


def bearing(la1, lo1, la2, lo2):
    p = np.pi / 180
    y = np.sin((lo2 - lo1) * p) * np.cos(la2 * p)
    x = np.cos(la1 * p) * np.sin(la2 * p) - np.sin(la1 * p) * np.cos(la2 * p) * np.cos((lo2 - lo1) * p)
    return (np.degrees(np.arctan2(y, x)) + 360) % 360


def nearest_place(lat, lon, keys=None):
    keys = keys or list(PLACES)
    d = {k: float(hav(lat, lon, PLACES[k][0], PLACES[k][1])) for k in keys}
    k = min(d, key=d.get)
    return k, d[k]


def operator_of(callsigns):
    for cs in callsigns:
        for pre, (name, code) in CALLSIGN_OPERATORS.items():
            if cs.startswith(pre) and cs[len(pre):].isdigit():
                return name, code
    return "", ""


def mins_after_1500(ts):
    """pandas Timestamp (MDT, naive) -> minutes after 15:00 MDT Oct. 17, 2026."""
    return (pd.Timestamp(ts) - pd.Timestamp("2026-10-17 15:00:00")).total_seconds() / 60.0


def clock(m):
    """minutes after 15:00 MDT Oct. 17 -> 'HH:MM' (24 h; '+1d' after midnight)."""
    if m is None or not np.isfinite(m):
        return ""
    t = pd.Timestamp("2026-10-17 15:00:00") + pd.Timedelta(minutes=round(float(m)))
    return t.strftime("%H:%M") + (" (+1d)" if t.day == 18 else "")


def front_time_at(lat, lon):
    """v5 front plane: minutes after 15:00 MDT at which the surface cold front passed (lat, lon). Gradient from front_plane.json."""
    fp = json.load(open(os.path.join(V5, "front_plane.json")))["fit"]
    from pyproj import Transformer
    tr = Transformer.from_crs(4326, 3400, always_xy=True)
    x0, y0 = tr.transform(FIRE[1], FIRE[0])
    x, y = tr.transform(np.asarray(lon), np.asarray(lat))
    b, c = fp["gradient_min_per_km"]
    return fp["t_fire_min"] + b * (x - x0) / 1000 + c * (y - y0) / 1000


def load_traces(icaos=None):
    """All positions from the part-1 adsb.lol trace files (both daily releases), one row per (icao, time).
    readsb trace point layout: [dt_s, lat, lon, alt_baro (ft or "ground"), gs_kt, track_deg, flags, vert_rate_fpm, details, source,
    geom_alt_ft, ...]."""
    rows = []
    for f in sorted(glob.glob(os.path.join(TRACES, "*", "trace_full_*.json"))):
        ic = os.path.basename(f)[len("trace_full_"):-5]
        if icaos is not None and ic not in icaos:
            continue
        raw = open(f, "rb").read()
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        tr = json.loads(raw)
        base, rel = tr.get("timestamp", 0), os.path.basename(os.path.dirname(f))
        mil = bool(tr.get("dbFlags", 0) & 1)
        for p in tr.get("trace", []):
            det = p[8] if len(p) > 8 and isinstance(p[8], dict) else {}
            rows.append((tr.get("icao"), tr.get("r", ""), tr.get("t", ""), tr.get("desc", ""), mil, base + p[0], p[1], p[2],
                         p[3], p[4], p[5], p[7] if len(p) > 7 else None, (det.get("flight") or "").strip(),
                         p[9] if len(p) > 9 else None, p[10] if len(p) > 10 else None, rel))
    A = pd.DataFrame(rows, columns=["icao", "reg", "type", "desc", "military_db", "t", "lat", "lon", "alt_baro", "gs_kt",
                                    "track_deg", "vrate_fpm", "callsign", "source", "geom_alt_ft", "release"])
    A = A.sort_values(["icao", "t"]).drop_duplicates(["icao", "t"]).reset_index(drop=True)
    A["mdt"] = pd.to_datetime(A.t, unit="s") - MDT_OFF
    A["on_ground"] = A.alt_baro.astype(str) == "ground"
    A["alt_ft"] = pd.to_numeric(A.alt_baro, errors="coerce")
    return A
