"""Polar-orbiter (VIIRS / MODIS) daytime passes over the CFB Wainwright fire, Oct. 17, 2026,
13:00-18:00 MDT (19:00-00:00 UTC), and VIIRS close-up imagery for each useful pass.

Usage:
  python3 -I polar_visible.py passes  OUT_DIR             # pass table (passes.csv), no network needed
  python3 -I polar_visible.py extract OUT_DIR [SAT ...]   # VIIRS SDR crops from NOAA Open Data on AWS
  python3 -I polar_visible.py render  OUT_DIR             # close-up PNGs from the crops

`passes` propagates the TLEs below (SGP4 via skyfield) and adds the scan angle, view zenith, pixel
size, sun elevation and the GOES-18 cloud indicators at the fire (from raw/goes/*_point_series.csv).

`extract` reads the VIIRS SDR HDF5 files lazily over anonymous HTTPS (fsspec + h5py) from
noaa-nesdis-{snpp,n20,n21}-pds: it reads the terrain-corrected GEO file (GITCO / GMTCO), finds the rows
and columns within ~110 km of the fire, then reads only that block of I1/I2/I3/I4 (and M3/M4/M5 for
true colour). No granule is written to disk. The crops are saved as crops_<sat>_<UTC>.npz.
STATUS (Oct. 27, 2026): `extract` and `render` have NOT been run. The three buckets were refused by
this session's egress proxy (HTTP 403 on CONNECT), so the code is untested against real granules.
Check bucket key layout and dataset names on first use.

Bands: I1 0.64 um, I2 0.865 um, I3 1.61 um (375 m), I4 3.74 um (375 m), M3 0.488, M4 0.555, M5 0.672 um
(750 m). Credit: NOAA/NASA VIIRS via NOAA Open Data on AWS.
"""
import csv
import datetime as dt
import os
import re
import sys

import numpy as np

FIRE = (52.7645, -110.9052)
TOWN = (52.8334, -110.8601)
T0, T1 = dt.datetime(2026, 10, 17, 19, 0), dt.datetime(2026, 10, 18, 0, 0)
MDT = dt.timedelta(hours=-6)
CREDIT = "NOAA/NASA VIIRS via NOAA Open Data on AWS"
R_E = 6371.0
GOES_RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "../../data/rubicon-ridge/raw/goes")

# TLEs fetched Oct. 27, 2026 from SatNOGS DB (db.satnogs.org, mirrors Space-Track). Checksums verified.
# Back-propagated 5-15 days to Oct. 17: expect pass-time errors of order seconds to ~1 min.
TLES = {
    "NOAA-20": ("1 43013U 17073A   26300.23278431  .00000035  00000-0  37611-4 0  9990",
                "2 43013  98.7892 238.5757 0000417 133.9273 226.1937 14.19533286463200"),
    "S-NPP":   ("1 37849U 11061A   26300.28937038  .00000041  00000-0  40343-4 0  9994",
                "2 37849  98.8082 240.3224 0001201 344.6360  15.4780 14.19541295777181"),
    "NOAA-21": ("1 54234U 22150A   26295.82296925  .00000030  00000-0  34999-4 0  9996",
                "2 54234  98.7164 232.5079 0001775 155.2159 204.9103 14.19554972204654"),
    "Aqua":    ("1 27424U 02022A   26276.25673623  .00000422  00000-0  93683-4 0  9993",
                "2 27424  98.4426 247.7760 0001135  62.4511  62.8039 14.62220671299061"),
    "Terra":   ("1 25994U 99068A   26286.28538608  .00000367  00000-0  83598-4 0  9998",
                "2 25994  97.9319 331.0709 0001867 180.2012 250.8881 14.61175622426850"),
}
# instrument, max scan angle (deg), nadir pixel (km) of the finest band
SENSOR = {"NOAA-20": ("VIIRS", 56.06, 0.375), "S-NPP": ("VIIRS", 56.06, 0.375),
          "NOAA-21": ("VIIRS", 56.06, 0.375), "Aqua": ("MODIS", 55.0, 0.25), "Terra": ("MODIS", 55.0, 0.25)}
BUCKET = {"S-NPP": ("noaa-nesdis-snpp-pds", "npp"), "NOAA-20": ("noaa-nesdis-n20-pds", "j01"),
          "NOAA-21": ("noaa-nesdis-n21-pds", "j02")}


def solar_elevation(t, lat, lon):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from goes_visible import solar_zenith
    return 90 - float(solar_zenith(t, lat, lon))


def pixel_km(sensor, scan_deg, h):
    """Approximate ground pixel size (along-scan, along-track) in km for the finest band.
    Along-track grows with slant range; along-scan with the ground-arc derivative ds/dtheta.
    VIIRS I-band: native along-scan sampling is aggregated 3:1 (scan < 31.59 deg), 2:1 (< 44.68 deg),
    1:1 beyond, which keeps pixels near 0.4-0.8 km. MODIS 250 m has no aggregation (bow-tie growth)."""
    th = np.radians(scan_deg)
    sin_g = (R_E + h) / R_E * np.sin(th)          # sin of view zenith at the ground
    vz = np.arcsin(np.clip(sin_g, -1, 1))
    slant = R_E * np.sin(vz - th) / np.sin(th) if th > 1e-6 else h
    dsdth = R_E * ((R_E + h) / R_E * np.cos(th) / np.cos(vz) - 1)   # km per rad
    if sensor == "VIIRS":
        agg = 3 if scan_deg < 31.59 else (2 if scan_deg < 44.68 else 1)
        native = 0.129 / h                         # rad per native along-scan sample (~0.129 km at nadir)
        return agg * native * dsdth, 0.371 * slant / h
    d = 0.25 / h
    return d * dsdth, 0.25 * slant / h


def goes_cloud(t):
    """GOES-18 10.3 um BT (fire pixel, burn-box median, background median) and FDC mask at the fire for
    the G18 Full Disk scan whose fire-latitude time (start + ~40 s) is closest to t."""
    try:
        import pandas as pd
        r = pd.read_csv(os.path.join(GOES_RAW, "rad_point_series.csv"), parse_dates=["scan_start_utc"])
        f = pd.read_csv(os.path.join(GOES_RAW, "fdc_point_series.csv"), parse_dates=["scan_start_utc"])
    except Exception:
        return {}
    m = r.merge(f[["scan_start_utc", "sat", "fire_Mask"]], on=["scan_start_utc", "sat"])
    m = m[m.sat == "G18"].copy()
    m["tdiff"] = (m.scan_start_utc + pd.Timedelta(seconds=40) - pd.Timestamp(t)).abs()
    row = m.sort_values("tdiff").iloc[0]
    if row["tdiff"] > pd.Timedelta(minutes=6):
        return {}
    return {"g18_scan_mdt": row.scan_start_mdt[11:16], "g18_c13_fire_k": round(row.c13_fire_pix_bt_k, 1),
            "g18_c13_burnbox_med_k": round(row.c13_burnbox_median_bt_k, 1),
            "g18_c13_bg_med_k": round(row.c13_bg_median_bt_k, 1), "g18_fdc_mask": int(row.fire_Mask)}


def passes(out):
    from skyfield.api import EarthSatellite, load, wgs84
    ts = load.timescale()
    site = wgs84.latlon(FIRE[0], FIRE[1], 700)
    rows = []
    for name, (l1, l2) in TLES.items():
        sat = EarthSatellite(l1, l2, name, ts)
        sec = np.arange(-3600, int((T1 - T0).total_seconds()) + 3600, 1.0)
        t = ts.utc(T0.year, T0.month, T0.day, T0.hour, 0, sec)
        alt, az, dist = (sat - site).at(t).altaz()
        el = alt.degrees
        for k in range(1, len(el) - 1):
            if not (el[k] > el[k - 1] and el[k] >= el[k + 1] and el[k] > 0):
                continue
            tu = t[k].utc_datetime().replace(tzinfo=None, microsecond=0)
            h = float(wgs84.height_of(sat.at(t[k])).km)
            sp = wgs84.subpoint_of(sat.at(t[k]))
            vz = 90 - el[k]
            scan = float(np.degrees(np.arcsin((R_E + 0.7) / (R_E + h) * np.sin(np.radians(vz)))))
            sensor, smax, _ = SENSOR[name]
            ps, pt = pixel_km(sensor, min(scan, smax), h)
            vel = sat.at(t[k + 1]).position.km - sat.at(t[k - 1]).position.km
            zpos = sat.at(t[k]).position.km
            node = "ascending" if vel[2] > 0 else "descending"
            sun = solar_elevation(tu, *FIRE)
            g = goes_cloud(tu)
            rows.append(dict(
                satellite=name, sensor=sensor, utc=tu.strftime("%Y-%m-%d %H:%M:%S"),
                mdt=(tu + MDT).strftime("%H:%M:%S"), in_window=T0 <= tu < T1, node=node,
                view_zenith_deg=round(vz, 1), scan_angle_deg=round(scan, 1), max_scan_deg=smax,
                in_swath=scan <= smax, sat_azimuth_from_fire_deg=round(float(az.degrees[k])),
                sat_alt_km=round(h), subsat_lat=round(float(sp.latitude.degrees), 2),
                subsat_lon=round(float(sp.longitude.degrees), 2),
                ground_dist_km=round(float(R_E * np.radians(vz - scan))),
                pix_along_scan_km=round(ps, 2) if scan <= smax else "", pix_along_track_km=round(pt, 2) if scan <= smax else "",
                sun_elev_deg=round(sun, 1), daytime=sun > 0, tle_age_days=round(float(t[k] - sat.epoch), 1), **g))
    rows.sort(key=lambda r: r["utc"])
    os.makedirs(out, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(out, "passes.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print({k: r.get(k) for k in ("satellite", "mdt", "scan_angle_deg", "in_swath", "sun_elev_deg",
                                     "pix_along_scan_km", "pix_along_track_km", "g18_c13_fire_k")})


# ---------------------------------------------------------------------------------------------- VIIRS
IBANDS = {"I1": "SVI01", "I2": "SVI02", "I3": "SVI03", "I4": "SVI04"}
MBANDS = {"M3": "SVM03", "M4": "SVM04", "M5": "SVM05"}
PREFIX = {"SVI01": "VIIRS-I1-SDR", "SVI02": "VIIRS-I2-SDR", "SVI03": "VIIRS-I3-SDR", "SVI04": "VIIRS-I4-SDR",
          "SVM03": "VIIRS-M3-SDR", "SVM04": "VIIRS-M4-SDR", "SVM05": "VIIRS-M5-SDR",
          "GITCO": "VIIRS-IMG-GEO-TC", "GMTCO": "VIIRS-MOD-GEO-TC"}
CROP_KM = 110


def granule_times(key):
    m = re.search(r"_d(\d{8})_t(\d{6})\d_e(\d{6})\d", key)
    d, a, b = m.groups()
    t0 = dt.datetime.strptime(d + a, "%Y%m%d%H%M%S")
    t1 = dt.datetime.strptime(d + b, "%Y%m%d%H%M%S")
    return t0, t1 + (dt.timedelta(days=1) if t1 < t0 else dt.timedelta(0))


def list_granules(bucket, prod, tag, t_pass):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import goes_common as gc
    keys = [k for k, _ in gc.list_keys(bucket, f"{PREFIX[prod]}/{t_pass:%Y/%m/%d}/")
            if os.path.basename(k).startswith(f"{prod}_{tag}_")]
    return [k for k in keys if granule_times(k)[0] - dt.timedelta(seconds=90) <= t_pass
            <= granule_times(k)[1] + dt.timedelta(seconds=90)]


def h5_open(fs, bucket, key):
    import h5py
    return h5py.File(fs.open(f"https://{bucket}.s3.amazonaws.com/{key}", "rb"), "r")


def read_scaled(f, band):
    """Return (scaled array reader, kind). SDR stores uint16 + [scale, offset] factors; fill >= 65528."""
    grp = [g for g in f["All_Data"] if g.endswith("SDR_All")][0]
    g = f["All_Data"][grp]
    var = "Reflectance" if "Reflectance" in g else "BrightnessTemperature"
    fac = g[var + "Factors"][:2] if var + "Factors" in g else np.array([1.0, 0.0])
    return g[var], fac, var


def extract(out, sats):
    import fsspec
    fs = fsspec.filesystem("https", block_size=2 ** 20)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "passes.csv")) as f:
        plist = [r for r in csv.DictReader(f) if r["sensor"] == "VIIRS" and r["in_window"] == "True"
                 and r["in_swath"] == "True" and r["daytime"] == "True" and r["satellite"] in sats]
    from pyproj import Geod
    geod = Geod(ellps="WGS84")
    for p in plist:
        bucket, tag = BUCKET[p["satellite"]]
        tp = dt.datetime.strptime(p["utc"], "%Y-%m-%d %H:%M:%S")
        for geo_prod, bands, res in (("GITCO", IBANDS, "I"), ("GMTCO", MBANDS, "M")):
            for gk in list_granules(bucket, geo_prod, tag, tp):
                with h5_open(fs, bucket, gk) as g:
                    gg = g["All_Data"][[x for x in g["All_Data"]][0]]
                    lat, lon = gg["Latitude"][:], gg["Longitude"][:]
                    ok = (lat > -90) & (lat < 90)
                    _, _, d = geod.inv(np.full(lat.shape, FIRE[1]), np.full(lat.shape, FIRE[0]),
                                       np.where(ok, lon, FIRE[1]), np.where(ok, lat, FIRE[0]))
                    near = ok & (d < CROP_KM * 1e3)
                    if not near.any():
                        continue
                    jj, ii = np.where(near)
                    J, I = slice(jj.min(), jj.max() + 1), slice(ii.min(), ii.max() + 1)
                    rec = dict(lat=lat[J, I], lon=lon[J, I], satz=gg["SatelliteZenithAngle"][J, I],
                               solz=gg["SolarZenithAngle"][J, I])
                t0g, _ = granule_times(gk)
                for bname, prod in bands.items():
                    bk = [k for k in list_granules(bucket, prod, tag, tp) if granule_times(k)[0] == t0g]
                    if not bk:
                        print("  missing", prod, t0g)
                        continue
                    with h5_open(fs, bucket, bk[0]) as b:
                        ds, fac, var = read_scaled(b, bname)
                        raw = ds[J, I]
                        val = raw.astype("float32") * fac[0] + fac[1]
                        val[raw >= 65528] = np.nan
                        rec[bname] = val
                fn = os.path.join(out, f"crops_{p['satellite']}_{res}_{t0g:%Y%m%dT%H%M%S}.npz")
                np.savez_compressed(fn, t_pass_utc=p["utc"], granule=os.path.basename(gk), **rec)
                print("wrote", fn, {k: v.shape for k, v in rec.items()}, flush=True)


# ---------------------------------------------------------------------------------------------- render
def render(out):
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pyproj import Proj
    from scipy.interpolate import griddata
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from goes_visible import stretch
    pr = Proj(proj="aeqd", lat_0=FIRE[0], lon_0=FIRE[1], ellps="WGS84")
    tx, ty = (v / 1e3 for v in pr(TOWN[1], TOWN[0]))
    for fn in sorted(glob.glob(os.path.join(out, "crops_*.npz"))):
        z = dict(np.load(fn))
        res = 0.375 if "_I_" in fn else 0.75
        x, y = pr(z["lon"], z["lat"])
        x, y = x / 1e3, y / 1e3
        sat = os.path.basename(fn).split("_")[1]
        tpass = dt.datetime.strptime(str(z["t_pass_utc"]), "%Y-%m-%d %H:%M:%S")
        satz = float(np.nanmedian(z["satz"][(np.hypot(x, y) < 5)])) if (np.hypot(x, y) < 5).any() else np.nan
        for half in (30, 75):
            g = np.arange(-half, half + res, res)
            GX, GY = np.meshgrid(g, g[::-1])
            sel = (np.abs(x) < half + 5) & (np.abs(y) < half + 5)
            pts = np.c_[x[sel], y[sel]]

            def grid(a, method="linear"):
                return griddata(pts, a[sel], (GX, GY), method=method)
            panels = []
            if "I1" in z:
                panels.append(("red", "I1 0.64 µm reflectance", stretch(grid(z["I1"]))))
            if all(b in z for b in ("M5", "M4", "M3")):
                panels.append(("truecolour", "true colour M5/M4/M3",
                               stretch(grid(z["M5"]), grid(z["M4"]), grid(z["M3"]))))
            hot = None
            if "I4" in z:
                from scipy.ndimage import median_filter
                bt = grid(z["I4"], "nearest")
                hot = (bt - median_filter(np.nan_to_num(bt), size=15)) > 10
                hot &= np.hypot(GX, GY) < 15
            for tag, label, img in panels:
                fig, ax = plt.subplots(figsize=(7, 7.4), dpi=110)
                ax.imshow(img, cmap="gray" if img.ndim == 2 else None, extent=(g[0], g[-1], g[0], g[-1]),
                          interpolation="nearest")
                if hot is not None and hot.any():
                    ax.contour(GX, GY, hot.astype(float), levels=[0.5], colors="red", linewidths=0.8)
                ax.plot(0, 0, "+", color="red", ms=14, mew=1.6)
                ax.plot(tx, ty, "s", mfc="none", mec="yellow", ms=7, mew=1.3)
                ax.annotate("Wainwright", (tx, ty), (tx + 1.5, ty + 1.5), color="yellow", fontsize=8)
                bx, by = g[0] + 0.06 * 2 * half, g[0] + 0.05 * 2 * half
                ax.plot([bx, bx + 10], [by, by], color="white", lw=3)
                ax.text(bx + 5, by + 0.02 * 2 * half, "10 km", color="white", ha="center", fontsize=9)
                ax.set_xlabel("km east of fire centre")
                ax.set_ylabel("km north of fire centre")
                ax.set_title(f"{sat} VIIRS {label}\n{tpass + MDT:%Y-%m-%d %H:%M} MDT ({tpass:%H:%M} UTC), "
                             f"view zenith {satz:.0f}°", fontsize=10)
                fig.text(0.02, 0.01, CREDIT + " | contrast enhanced per image" +
                         (" | red outline: I4 hot pixels" if hot is not None and hot.any() else ""), fontsize=7)
                fig.tight_layout(rect=(0, 0.02, 1, 1))
                pf = os.path.join(out, f"{sat}_{tag}_{2 * half}km_{tpass:%Y%m%dT%H%M}Z.png")
                fig.savefig(pf)
                plt.close(fig)
                print("wrote", pf)


if __name__ == "__main__":
    cmd, out = sys.argv[1], sys.argv[2]
    if cmd == "passes":
        passes(out)
    elif cmd == "extract":
        extract(out, sys.argv[3:] or list(BUCKET))
    elif cmd == "render":
        render(out)
