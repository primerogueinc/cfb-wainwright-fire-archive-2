"""Rubicon Ridge v5, step 1: the moving front and the post-front wind / spread schedules.
Usage: python3 -I v5_inputs.py            (reads data/rubicon-ridge/raw; writes data/rubicon-ridge/model_v5/{front_plane,wind_schedules}.json
                                           and schedules.csv)

1. Front plane. Arrival time of the cold front t = a + b*x + c*y (x east, y north, km from the fire centre) fitted by least squares to
   the station / sensor times listed in data/rubicon-ridge/raw/README_pull.md and raw/acis/README.md (minutes after 15:00 MDT).
   The gradient (b, c) gives the front's direction of travel and speed; a is the time at the fire. Leave-one-out fits give the
   spread. The model keeps the fitted gradient and shifts a to 15:27, 15:29 and 15:31 (the measured window at the fire).
2. Post-front wind (10 m SUSTAINED means, never gusts). Five schedules, each a function of minutes since the front (tau):
   W1_czwn_obs     CZWN as observed: direction 290 -> 300 -> 307 -> 304 -> 299, speed 61 km/h for the first half hour (the
                   15:16-16:00 vector mean after removing the 16 pre-front minutes), 47.6 at the 16:00 10-min mean, then the
                   hourly means 45.2 / 37.2 / 30.3 / 24.8 km/h. Times are relative to the front at CZWN (15:16) and the pattern is
                   carried with the front to the fire.
   W2_fixed300     300 deg, 48 km/h (CZWN 10-min mean at 16:00, the "45-50 km/h") to tau=74 min, then the CZWN hourly decline.
   W3_fixed310     as W2 at 310 deg (Hughenden, the nearest station, 308-312 deg after the front).
   W4_fixed315     as W2 at 315 deg (v4's best-fitting direction).
   W5_models       RAP (13 km) and HRDPS (2.5 km, 12Z and 18Z runs) 10 m wind at the fire: model mean of RAP and the HRDPS mean,
                   at 16:00, 17:00, 18:00 and 19:00 MDT (absolute clock times).
3. Fine-fuel moisture after the front: hourly FFMC (Van Wagner 1977) run on CZWN hourly climate data from CWFIS daily FFMC 87.2
   at 1 p.m. Oct. 16 (as part 8b / v4), continued through the evening (cooler, 40-50% humidity after the front).
"""
import os, sys, json, math
import numpy as np, pandas as pd
from pyproj import Transformer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v5_common import REPO, V4REPO, OUT, FIRE_LONLAT, hffmc, isi, ros_grass, lb_grass
os.makedirs(OUT, exist_ok=True)
RAW = os.path.join(REPO, "data/rubicon-ridge/raw")

# ---------------------------------------------------------------- 1. front plane
ST = [  # name, lat, lon, minutes after 15:00 MDT, source
    ("CFB Wainwright CZWN", 52.8294, -111.0956, 16, "ECCC hourly SWOB, peak gust 15:16 (82 km/h)"),
    ("Albert Hall AGCM", 53.1134, -111.1812, 2, "ACIS hourly humidity (raw/acis/README.md), ~15:00-15:05"),
    ("Bellshill AGCM", 52.5824, -111.4654, 2, "ACIS hourly humidity, ~15:00-15:05"),
    ("Gilt Edge North AGCM", 53.0710, -110.6210, 30, "ACIS hourly humidity ~15:30 (peak gust 15:28)"),
    ("Hughenden AGCM", 52.5782, -110.7837, 39, "ACIS hourly humidity ~15:39 (peak-gust proxy 15:43)"),
    ("Wainwright town (PurpleAir)", 52.8334, -110.8601, 28, "PurpleAir dust spike 15:28"),
    ("Castor (PurpleAir)", 52.2197, -111.9063, 7, "PurpleAir dust spike 15:07"),
    ("Coronation CWCT", 52.07427, -111.449205, 31, "ECCC 1-min, wind shift 15:31"),
    ("Lloydminster CYLL", 53.3106, -110.0740, 45, "ECCC SPECIs 15:41-15:48, ~15:45"),
]
tf = Transformer.from_crs(4326, 3400, always_xy=True).transform
fx, fy = tf(*FIRE_LONLAT)
X = np.array([[(tf(lo, la)[0] - fx) / 1000, (tf(lo, la)[1] - fy) / 1000] for _, la, lo, _, _ in ST]); t = np.array([s[3] for s in ST], float)
def fit(idx):
    A = np.c_[np.ones(len(idx)), X[idx]]; c, *_ = np.linalg.lstsq(A, t[idx], rcond=None); return c
def desc(c):
    g = c[1:]; return dict(t_fire_min=round(float(c[0]), 2), speed_kmh=round(60 / float(np.hypot(*g)), 1),
                           toward_deg=round(float(np.degrees(np.arctan2(g[0], g[1])) % 360), 1))
call = fit(np.arange(len(ST))); res = t - np.c_[np.ones(len(ST)), X] @ call
loo = [desc(fit(np.delete(np.arange(len(ST)), k))) | {"left_out": ST[k][0]} for k in range(len(ST))]
front = dict(method="least-squares plane t = a + b*x + c*y (km, EPSG:3400, origin at the fire centre)", fire_centre_lonlat=FIRE_LONLAT,
             stations=[dict(name=s[0], lat=s[1], lon=s[2], t_mdt=f"15:{s[3]:02d}", source=s[4], residual_min=round(float(r), 1)) for s, r in zip(ST, res)],
             fit=desc(call) | dict(gradient_min_per_km=[round(float(call[1]), 4), round(float(call[2]), 4)], rms_residual_min=round(float(np.sqrt(np.mean(res ** 2))), 1)),
             leave_one_out=loo,
             loo_range=dict(t_fire=[min(d["t_fire_min"] for d in loo), max(d["t_fire_min"] for d in loo)],
                            toward=[min(d["toward_deg"] for d in loo), max(d["toward_deg"] for d in loo)],
                            speed=[min(d["speed_kmh"] for d in loo), max(d["speed_kmh"] for d in loo)]),
             front_times_at_fire_used=["15:27", "15:29", "15:31"],
             note="Model uses the fitted gradient (direction and speed) and shifts the time at the fire to 15:27 / 15:29 / 15:31.")
json.dump(front, open(os.path.join(OUT, "front_plane.json"), "w"), indent=1)
print(json.dumps(front["fit"]), front["loo_range"], np.round(res, 1))

# ---------------------------------------------------------------- 3. FFMC (hourly)
c = pd.read_csv(os.path.join(V4REPO, "data/weather_exports/3_climate_hourly_czwn_oct2026.csv"))
c["t"] = (pd.to_datetime(c["time_mdt_local_clock_utc-6"]) - pd.Timestamp("2026-10-17")).dt.total_seconds() / 3600
c = c[(c.t >= -11) & (c.t <= 21)].sort_values("t")
ff = 87.2; FF = []
for _, r in c.iterrows():
    if r.t >= -11 + 1e-6: ff = hffmc(ff, r["Temp (°C)"], r["Rel Hum (%)"], r["Wind Spd (km/h)"])
    FF.append((r.t, ff))
FF = np.array(FF)
def ffmc_at(tmin):   # minutes after 15:00 MDT -> FFMC (hourly values, linear between)
    return float(np.interp(15 + tmin / 60, FF[:, 0], FF[:, 1]))

# ---------------------------------------------------------------- 2. wind schedules
sw = pd.read_csv(os.path.join(RAW, "eccc_swob/swob_hourly_key_vars.csv"))
cz = sw[(sw.station == "CZWN")].set_index("obs_time_mdt")
def vec(s, d): a = math.radians(d); return np.array([s * math.sin(a), s * math.cos(a)])
def unvec(v): return float(np.hypot(*v)), float(math.degrees(math.atan2(v[0], v[1])) % 360)
h16 = cz.loc["2026-10-17 16:00"]; h15 = cz.loc["2026-10-17 15:00"]
pre = vec(h15.avg_wnd_spd_10m_pst10mts, h15.avg_wnd_dir_10m_pst10mts)          # pre-front wind, 14:50-15:00 mean
first = unvec((60 * vec(h16.avg_wnd_spd_10m_pst1hr, h16.avg_wnd_dir_10m_pst1hr) - 16 * pre) / 44)   # 15:16-16:00 mean
hourly = {hh: (cz.loc[f"2026-10-17 {hh}:00"].avg_wnd_spd_10m_pst1hr, cz.loc[f"2026-10-17 {hh}:00"].avg_wnd_dir_10m_pst1hr) for hh in (17, 18, 19, 20)}
m10 = (h16.avg_wnd_spd_10m_pst10mts, h16.avg_wnd_dir_10m_pst10mts)
# tau relative to the front at CZWN (15:16): 16:00 = 44, 16:30 = 74 (17:00 hourly mean), 17:30 = 134, 18:30 = 194, 19:30 = 254
cz_spd = [(0, first[0]), (30, first[0]), (44, m10[0]), (74, hourly[17][0]), (134, hourly[18][0]), (194, hourly[19][0]), (254, hourly[20][0])]
cz_dir = [(0, 290.0), (10, first[1]), (44, m10[1]), (74, hourly[17][1]), (134, hourly[18][1]), (194, hourly[19][1]), (254, hourly[20][1])]
flat_spd = [(0, 48.0), (74, hourly[17][0]), (134, hourly[18][0]), (194, hourly[19][0]), (254, hourly[20][0])]
rap = pd.read_csv(os.path.join(RAW, "rap/rap_points_f00.csv")); rap = rap[rap.point == "fire_centre"].set_index("valid_mdt")
hr = pd.read_csv(os.path.join(RAW, "eccc_models/hrdps_points_wide.csv")); hr = hr[(hr.point == "fire_centre") & hr.run.isin(["12Z", "18Z"])]
mod = []
for hh in (16, 17, 18, 19):
    key = f"2026-10-17 {hh}:00"
    vr = vec(rap.loc[key, "wspd10m_kmh_idw"], rap.loc[key, "wdir10m_deg_idw"])
    hh_rows = hr[hr.valid_mdt == key]
    vh = np.mean([vec(r["WIND_AGL-10m"] * 3.6, r["WDIR_AGL-10m"]) for _, r in hh_rows.iterrows()], axis=0) if len(hh_rows) else vr
    s_, d_ = unvec((vr + vh) / 2); mod.append((hh, s_, d_, unvec(vr), unvec(vh)))
SCHED = {
    "W1_czwn_obs": dict(kind="tau_station", spd=cz_spd, dir=cz_dir, note="CZWN observed (veer 290->300->307, 61 -> 48 -> 45 -> 37 km/h)"),
    "W2_fixed300": dict(kind="tau_station", spd=flat_spd, dir=[(0, 300.0), (254, 300.0)], note="300 deg, 48 km/h then CZWN hourly decline"),
    "W3_fixed310": dict(kind="tau_station", spd=flat_spd, dir=[(0, 310.0), (254, 310.0)], note="310 deg (Hughenden 308-312), 48 km/h then decline"),
    "W4_fixed315": dict(kind="tau_station", spd=flat_spd, dir=[(0, 315.0), (254, 315.0)], note="315 deg (v4 best fit), 48 km/h then decline"),
    "W5_models": dict(kind="clock", spd=[((h - 15) * 60, s) for h, s, d, _, _ in mod], dir=[((h - 15) * 60, d) for h, s, d, _, _ in mod],
                      note="mean of RAP and HRDPS (12Z/18Z) 10 m wind at the fire, 16:00-19:00 MDT"),
}
def sample(name, t_fire, fuel, mult, shape, step=0.25, n_min=300):
    """Arrays over tau = 0..n_min (minutes since the front reached the fire): head direction (rad, direction the fire is pushed),
    head ROS (m/min), eccentricity; plus the wind speed, direction, FFMC, ISI used."""
    s = SCHED[name]; tau = np.arange(0, n_min + 1e-9, step)
    def interp(knots, x):
        k = np.array(knots, float); xx = np.asarray(x, float)
        if s["kind"] == "clock": xx = np.maximum(xx + t_fire, k[0, 0])       # absolute clock (minutes after 15:00)
        return np.interp(xx, k[:, 0], k[:, 1])
    ws = interp(s["spd"], tau)
    d = s["dir"]; a = np.radians(np.array(d, float)[:, 1])
    ux = interp([(k[0], math.sin(v)) for k, v in zip(d, a)], tau); uy = interp([(k[0], math.cos(v)) for k, v in zip(d, a)], tau)
    wd = np.degrees(np.arctan2(ux, uy)) % 360
    ffm = np.array([ffmc_at(t_fire + x) for x in tau]); I = np.array([isi(f, w) for f, w in zip(ffm, ws)])
    ros = mult * np.array([ros_grass(i, fuel) for i in I])
    lb = np.array([lb_grass(w) for w in ws]) if shape == "steady" else np.full(len(tau), 3.0)
    ecc = np.sqrt(1 - 1 / lb ** 2)
    head = np.radians((wd + 180) % 360)
    return dict(tau=tau, ws=ws, wd=wd, ffmc=ffm, isi=I, ros=ros, ecc=ecc, head=head)

if __name__ == "__main__":
    json.dump(dict(schedules={k: {kk: vv for kk, vv in v.items()} for k, v in SCHED.items()},
                   czwn_first_post_front_segment=dict(speed_kmh=round(first[0], 1), dir_deg=round(first[1]), method="(60 min hourly mean vector - 16 min x 14:50-15:00 10-min mean) / 44"),
                   model_at_fire=[dict(mdt=f"{h}:00", mean_kmh=round(s_, 1), mean_dir=round(d_), rap=(round(r[0], 1), round(r[1])), hrdps=(round(q[0], 1), round(q[1]))) for h, s_, d_, r, q in mod],
                   ffmc_hourly=[dict(hour_mdt=f"{int(x) % 24:02d}:00", ffmc=round(y, 1)) for x, y in FF if x >= 12],
                   curing_pct=95, note="Sustained 10 m means only; gusts are not used for spread."),
              open(os.path.join(OUT, "wind_schedules.json"), "w"), indent=1, default=float)
    rows = []
    for name in SCHED:
        for fuel in ("O-1a", "O-1b"):
            z = sample(name, 29.0, fuel, 1.0, "steady")
            for tt in (0, 10, 30, 60, 90, 120, 180):
                k = int(tt / 0.25)
                rows.append(dict(schedule=name, fuel=fuel, minutes_after_front=tt, wind_kmh=round(z["ws"][k], 1), wind_from_deg=round(z["wd"][k]),
                                 ffmc=round(z["ffmc"][k], 1), isi=round(z["isi"][k], 1), head_ros_kmh=round(z["ros"][k] * 0.06, 2)))
    df = pd.DataFrame(rows); df.to_csv(os.path.join(OUT, "schedules.csv"), index=False)
    print(df[df.minutes_after_front.isin([0, 30, 60, 120])].to_string(index=False))
