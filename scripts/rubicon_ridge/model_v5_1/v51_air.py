"""Rubicon Cloud (working model v5.1), step 1: the aircraft timeline and its link to model v5.
Usage:  python3 v51_air.py          (needs the part-1 checkout; set CFB_ARCHIVE_1 if it is not next to this repository)
Writes data/rubicon-ridge/model_v5_1_rubicon_cloud/:
  provenance_check.json         the supplied file vs the adsb.lol trace files in part 1 (row match, "base" point, filter)
  aircraft_summary.csv          one row per aircraft in the supplied file: first/last, closest approaches, call signs, operator
  flight_pieces.csv             every continuously-heard piece of each aircraft's flights, 12:00 MDT Oct. 17 - 06:00 Oct. 18,
                                with where it starts and ends (inferred)
  wainwright_visits.csv         arrivals into / departures from the Wainwright area, with estimated landing / take-off times
  response_timeline.csv         chronological events, with minutes after the front at the fire and after the v5 arrival medians
  airspace_1400_1800.csv        every aircraft heard within 80 km of the fire 14:00-18:00 MDT (any altitude), for the null check
  context_tracks.csv            full positions of the 11 aircraft, 12:00 MDT Oct. 17 - 06:00 Oct. 18 (from the part-1 archive)
  proximity_burn_downwind.csv   positions within 25 km of the fire: distance to the burn outline, bearing, downwind sector flag
  approach_paths.csv            per Wainwright arrival/departure: closest point to the burn, cross-track deviation, downwind entry
  vertical_irregularity.csv     altitude reversals and vertical-rate scatter near Wainwright vs near other airfields
  front_passage_c510.csv        the one aircraft aloft near the fire at front passage: its position relative to the v5 front plane
  air_timing_check.json         the air-informed timing bound and likelihood weights over the v5 runs; updated arrival summary
  run_weights_air.csv           v5 run weights x air likelihood (default scar and Satellogic refit)
"""
import os, sys, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import v51_common as C

os.makedirs(C.OUT, exist_ok=True)
J = {}                                     # everything that goes into the JSON summaries
pd.options.mode.chained_assignment = None

# ---------------------------------------------------------------- 1. provenance of the supplied file
raw_bytes = open(C.SUPPLIED, "rb").read()
S = pd.read_csv(C.SUPPLIED)
S["ts"] = (pd.to_datetime(S.utc) - pd.Timestamp("1970-01-01")).dt.total_seconds().round().astype("int64")
S["mdt_ts"] = pd.to_datetime(S.mdt)
A = C.load_traces()                        # all 151 aircraft in the part-1 pull
hexes = sorted(S.hex.unique())
best, err, alt_eq = [], [], []
for i, r in S.iterrows():
    c = A[(A.icao == r.hex) & ((A.t - r.ts).abs() <= 1.0)]
    d = (c.lat - r.lat).abs() + (c.lon - r.lon).abs()
    j = d.idxmin()
    best.append(j); err.append(max(abs(A.at[j, "lat"] - r.lat), abs(A.at[j, "lon"] - r.lon)))
    alt_eq.append(str(A.at[j, "alt_baro"]) == str(r.alt_ft))
S["archive_idx"] = best
# what does km_from_base measure? least-squares fit of a point to the given distances
from scipy.optimize import least_squares
fit = least_squares(lambda x: C.hav(S.lat.values, S.lon.values, x[0], x[1]) - S.km_from_base.values, [52.78, -110.9])
base = (round(float(fit.x[0]), 4), round(float(fit.x[1]), 4))
resid = C.hav(S.lat.values, S.lon.values, 52.80, -110.90) - S.km_from_base.values
# reconstruct the selection: within 40 km of the base point, numeric baro altitude <= 12,000 ft, in a time window
A["km_base"] = C.hav(A.lat.values, A.lon.values, 52.80, -110.90)
sel = (A.km_base <= 40.05) & (A.alt_ft <= 12000)
t0, t1 = S.ts.min(), S.ts.max()
inwin = sel & (A.t >= t0 - 0.5) & (A.t <= t1 + 0.5)
prev_q = A[sel & (A.t < t0 - 0.5)].t.max(); next_q = A[sel & (A.t > t1 + 0.5)].t.max() if (sel & (A.t > t1 + 0.5)).any() else None
next_q = A[sel & (A.t > t1 + 0.5)].t.min() if (sel & (A.t > t1 + 0.5)).any() else None
ts2mdt = lambda t: (pd.Timestamp(t, unit="s") - C.MDT_OFF).strftime("%Y-%m-%d %H:%M:%S")
in_s = set(S.archive_idx)
J["provenance"] = {
    "file": os.path.relpath(C.SUPPLIED, C.REPO), "sha256": hashlib.sha256(raw_bytes).hexdigest(), "bytes": len(raw_bytes),
    "rows": int(len(S)), "aircraft": int(S.hex.nunique()), "types": sorted(S.type.unique().tolist()),
    "columns": list(pd.read_csv(C.SUPPLIED, nrows=0).columns),
    "mdt_column_is_utc_minus_6h": bool(((pd.to_datetime(S.utc) - pd.to_datetime(S.mdt)) == pd.Timedelta(hours=6)).all()),
    "first_mdt": S.mdt.min(), "last_mdt": S.mdt.max(),
    "alt_ft_range": [int(S.alt_ft.min()), int(S.alt_ft.max())], "km_from_base_range": [float(S.km_from_base.min()), float(S.km_from_base.max())],
    "military_flag_values": sorted(S.military.astype(str).unique().tolist()), "on_ground_values": sorted(S.on_ground.astype(str).unique().tolist()),
    "match_to_adsblol_traces": {
        "archive": "part 1 data/adsb/oct17/raw_traces/ (adsb.lol globe_history_2026 releases v2026.10.17 and v2026.10.18, readsb prod)",
        "rows_matched_within_1s_same_icao": int(len(S)),
        "distinct_archive_points_matched": int(len(in_s)),
        "duplicate_hex_utc_rows_in_file": int(S.duplicated(["hex", "utc"]).sum()),
        "max_position_difference_deg": round(float(max(err)), 5),
        "baro_altitude_identical_share": round(float(np.mean(alt_eq)), 4),
        "note": "Times in the file are the archive's fractional-second times rounded to whole seconds; where two archive points round "
                "to the same second the nearest one is used, which explains the small position differences."},
    "km_from_base": {"base_point_fitted_latlon": base, "base_point_reading": "52.80 N, 110.90 W: a round-number point 4.0 km north of the v5 fire centre and "
                     "2.8 km west of Wainwright Airport (CYWV); not the CFB Wainwright airfield / CZWN (13.6 km away)",
                     "max_abs_residual_km_vs_52.80_-110.90": round(float(np.abs(resid).max()), 3)},
    "reconstructed_selection": {
        "rule": "positions within 40 km of 52.80 N 110.90 W, numeric barometric altitude <= 12,000 ft (excludes 'ground'), "
                "between a start time in [prev qualifying point, first row] and an end time in [last row, next qualifying point]",
        "archive_points_meeting_rule_in_window": int(inwin.sum()),
        "of_which_in_file": int((inwin & A.index.isin(in_s)).sum()),
        "file_rows_outside_rule": int((~A.loc[list(in_s)].index.isin(A[inwin].index)).sum()),
        "start_bracket_mdt": [ts2mdt(prev_q) if pd.notna(prev_q) else None, S.mdt.min()],
        "end_bracket_mdt": [S.mdt.max(), ts2mdt(next_q) if next_q else None],
        "reading": "Every qualifying archive point between the first and last row is in the file, so the file is this filter applied "
                   "to the adsb.lol archive. The start cut lies somewhere between the previous qualifying point (11:20:36 MDT) and the "
                   "first row (15:31:28 MDT), consistent with a cut at about 15:30 MDT; the end cut lies between 03:57:58 and 11:49:36 "
                   "MDT Oct. 18. Nothing below 4,250 ft appears because no position below that was received near the base, not "
                   "because of a floor in the filter."},
    "part1_already_has_these_positions": True,
    "provenance_statement": "Supplied by Prime Rogue; every row matches the adsb.lol open archive (ODbL 1.0) trace files already "
                            "archived in part 1 (data/adsb/oct17/raw_traces/). Data: (c) adsb.lol contributors, ODbL 1.0.",
}
print(json.dumps(J["provenance"], indent=1)[:2500])

# ---------------------------------------------------------------- 2. context tracks (the 11 aircraft) and airspace null check
W0, W1 = pd.Timestamp("2026-10-17 12:00"), pd.Timestamp("2026-10-18 06:00")
T = A[A.icao.isin(hexes) & (A.mdt >= W0) & (A.mdt <= W1)].copy()
T["km_fire"] = C.hav(T.lat.values, T.lon.values, *C.FIRE)
T["km_cywv"] = C.hav(T.lat.values, T.lon.values, *C.PLACES["CYWV Wainwright Airport"][:2])
T["in_supplied_file"] = T.index.isin(in_s)
cs_by = T[T.callsign != ""].groupby("icao").callsign.agg(lambda s: sorted(set(s)))
ctx = T[["icao", "reg", "type", "callsign", "mdt", "lat", "lon", "alt_baro", "gs_kt", "track_deg", "vrate_fpm", "geom_alt_ft",
         "km_fire", "km_cywv", "in_supplied_file", "release"]].copy()
ctx["mdt"] = ctx.mdt.dt.strftime("%Y-%m-%d %H:%M:%S.%f").str[:-4]
ctx["lat"] = ctx.lat.round(5); ctx["lon"] = ctx.lon.round(5); ctx["km_fire"] = ctx.km_fire.round(2); ctx["km_cywv"] = ctx.km_cywv.round(2)
ctx.to_csv(os.path.join(C.OUT, "context_tracks.csv"), index=False)

B = A[(A.mdt >= "2026-10-17 14:00") & (A.mdt <= "2026-10-17 18:00")].copy()
B["km_fire"] = C.hav(B.lat.values, B.lon.values, *C.FIRE)
B = B[B.km_fire <= 80]
air = B.groupby(["icao", "reg", "type"]).agg(callsigns=("callsign", lambda s: " ".join(sorted(set(s) - {""}))),
                                              military_db=("military_db", "max"), first_mdt=("mdt", "min"), last_mdt=("mdt", "max"),
                                              min_km_fire=("km_fire", "min"), min_alt_ft=("alt_ft", "min"), max_alt_ft=("alt_ft", "max"),
                                              n=("t", "size")).reset_index()
air["below_12000_ft_within_80km"] = air.min_alt_ft <= 12000
air["first_mdt"] = air.first_mdt.dt.strftime("%H:%M:%S"); air["last_mdt"] = air.last_mdt.dt.strftime("%H:%M:%S")
air = air.sort_values("first_mdt"); air["min_km_fire"] = air.min_km_fire.round(1)
air.to_csv(os.path.join(C.OUT, "airspace_1400_1800.csv"), index=False)
low = B[B.alt_ft <= 12000]
J["airspace_null_check"] = {
    "window": "14:00-18:00 MDT, within 80 km of the fire centre, all aircraft in the part-1 pull",
    "aircraft_heard": int(len(air)), "aircraft_below_12000ft": int(air.below_12000_ft_within_80km.sum()),
    "military_flagged_aircraft_anywhere_in_pull": int(A.groupby("icao").military_db.max().sum()),
    "below_12000ft_by_aircraft": {r.reg: [r.first_mdt, r.last_mdt] for r in air[air.below_12000_ft_within_80km].itertuples()},
    "gap_below_12000ft_mdt": None,
}
lo = low.sort_values("t")
gaps = lo.mdt.diff().dt.total_seconds() / 60
k = gaps.idxmax() if len(lo) > 1 else None
if k is not None:
    J["airspace_null_check"]["gap_below_12000ft_mdt"] = [lo.mdt.shift(1)[k].strftime("%H:%M:%S"), lo.mdt[k].strftime("%H:%M:%S"),
                                                         round(float(gaps[k]), 1)]
print(json.dumps(J["airspace_null_check"], indent=1))

# ---------------------------------------------------------------- 3. flight pieces (continuously heard) -> flights, with start/end inference
def trend(q, first=True, minutes=3):
    w = q[q.mdt <= q.mdt.iloc[0] + pd.Timedelta(minutes=minutes)] if first else q[q.mdt >= q.mdt.iloc[-1] - pd.Timedelta(minutes=minutes)]
    a = w.alt_ft.dropna()
    if len(a) < 2:
        return 0.0, np.nan
    dt = (w.mdt.iloc[-1] - w.mdt.iloc[0]).total_seconds() / 60
    return float(a.iloc[-1] - a.iloc[0]), (float(a.iloc[-1] - a.iloc[0]) / dt if dt > 0.2 else np.nan)

HELI_TYPES = {"EC45", "BK17", "EC35", "B06", "B407", "AS50", "A139", "S76"}
def keys_for(actype):
    """Fixed-wing aircraft can only be landing at / leaving an aerodrome; helicopters may use a town or hospital helipad."""
    return list(C.PLACES) if actype in HELI_TYPES else [k for k in C.PLACES if k != "Wainwright town" and "Hospital" not in k
                                                         and "Medical" not in k]

def toward(lat, lon, trk, outbound=False, maxdeg=20.0, maxkm=60.0, keys=None):
    """The listed place the aircraft is heading toward (inbound) or appears to have come from (outbound): within maxkm, with the
    bearing within maxdeg of the ground track; the one with the smallest cross-track offset wins."""
    best = None
    for k in (keys or list(C.PLACES)):
        pl, pn, _ = C.PLACES[k]
        d = float(C.hav(lat, lon, pl, pn))
        if d > maxkm:
            continue
        b = float(C.bearing(pl, pn, lat, lon)) if outbound else float(C.bearing(lat, lon, pl, pn))
        dev = abs((b - trk + 180) % 360 - 180)
        if d < 3 or dev <= maxdeg:
            score = d * np.sin(np.radians(dev)) + 0.1 * d
            if best is None or score < best[3]:
                best = (k, d, dev, score)
    return best

def edge(q, first=True):
    """First (or last) point, with the track and ground speed taken from the nearest point within 60 s that has them."""
    r = (q.iloc[0] if first else q.iloc[-1]).copy()
    w = q[(q.mdt - r.mdt).abs() <= pd.Timedelta(seconds=60)].dropna(subset=["track_deg", "gs_kt"])
    if len(w):
        rr = w.iloc[0] if first else w.iloc[-1]
        r["track_deg"], r["gs_kt"] = rr.track_deg, rr.gs_kt
    return r

def classify_start(q, ground_before):
    KEYS = keys_for(q.type.iloc[0])
    r = edge(q, True); k, d = C.nearest_place(r.lat, r.lon, KEYS); elev = C.PLACES[k][2]
    climb, rate = trend(q, True)
    if ground_before or (r.alt_ft - elev <= 2000 and d <= 15):
        return "take-off (heard)", k, d, r.mdt
    def dev_from(kk, outbound):
        if not np.isfinite(r.track_deg):
            return 0.0
        pl, pn, _ = C.PLACES[kk]
        b = float(C.bearing(pl, pn, r.lat, r.lon)) if outbound else float(C.bearing(r.lat, r.lon, pl, pn))
        return abs((b - r.track_deg + 180) % 360 - 180)
    if d <= 8 and r.alt_ft - elev <= 4500 and climb >= -100 and (r.gs_kt < 100 or dev_from(k, True) <= 90):
        gs = r.gs_kt if r.gs_kt > 40 else 120.0
        m = max(max(r.alt_ft - elev, 0) / 1500.0, d / (gs * 1.852) * 60)
        return "first heard low near place (take-off estimated)", k, d, r.mdt - pd.Timedelta(minutes=m)
    if climb >= 300 and np.isfinite(r.track_deg):
        tw = toward(r.lat, r.lon, r.track_deg, outbound=True, maxdeg=30, keys=KEYS)
        if tw and r.alt_ft - C.PLACES[tw[0]][2] <= 12000:
            rr = max(rate if np.isfinite(rate) else 0.0, 500.0); gs = r.gs_kt if r.gs_kt > 60 else 150.0
            m = max((r.alt_ft - C.PLACES[tw[0]][2]) / rr, tw[1] / (gs * 1.852) * 60)
            return "climbing out (take-off estimated)", tw[0], tw[1], r.mdt - pd.Timedelta(minutes=m)
    if d <= 35 and r.alt_ft - elev <= 4500 and climb >= -100 and np.isfinite(r.track_deg):
        tw = toward(r.lat, r.lon, r.track_deg, outbound=True, maxdeg=45, maxkm=35, keys=KEYS)
        if tw:
            gs = r.gs_kt if r.gs_kt > 40 else 120.0
            m = max(max(r.alt_ft - C.PLACES[tw[0]][2], 0) / 1000.0, tw[1] / (gs * 1.852) * 60)
            return "first heard low, heading away (take-off estimated)", tw[0], tw[1], r.mdt - pd.Timedelta(minutes=m)
    return "first heard en route", k, d, pd.NaT

def classify_end(q, ground_after):
    KEYS = keys_for(q.type.iloc[0])
    r = edge(q, False); k, d = C.nearest_place(r.lat, r.lon, KEYS); elev = C.PLACES[k][2]
    ch, rate = trend(q, False)
    if ground_after or (r.alt_ft - elev <= 2000 and d <= 15):
        return "landing (heard)", k, d, r.mdt
    def dev_to(kk):
        if not np.isfinite(r.track_deg):
            return 0.0
        pl, pn, _ = C.PLACES[kk]
        return abs((float(C.bearing(r.lat, r.lon, pl, pn)) - r.track_deg + 180) % 360 - 180)
    if d <= 8 and r.alt_ft - elev <= 3500 and ch <= 100 and (r.gs_kt < 100 or dev_to(k) <= 60):
        gs = r.gs_kt if r.gs_kt > 40 else 120.0
        m = max(max(r.alt_ft - elev, 0) / 700.0, d / (gs * 1.852) * 60)
        return "last heard low near place (landing estimated)", k, d, r.mdt + pd.Timedelta(minutes=m)
    if (ch <= -200 or r.alt_ft - elev <= 4500) and np.isfinite(r.track_deg):
        tw = toward(r.lat, r.lon, r.track_deg, keys=KEYS)
        if tw:
            rr = max(-rate if np.isfinite(rate) else 0.0, 500.0); gs = r.gs_kt if r.gs_kt > 40 else 120.0
            m = max((r.alt_ft - C.PLACES[tw[0]][2]) / rr, tw[1] / (gs * 1.852) * 60)
            return "descending, lost below coverage (landing estimated)", tw[0], tw[1], r.mdt + pd.Timedelta(minutes=m)
    if d <= 35 and r.alt_ft - elev <= 4500 and ch <= 100 and np.isfinite(r.track_deg):
        tw = toward(r.lat, r.lon, r.track_deg, maxdeg=60, maxkm=35, keys=KEYS) or ((k, d, np.nan, 0.0) if ch <= -300 else None)
        if tw:
            gs = r.gs_kt if r.gs_kt > 40 else 120.0
            m = max(max(r.alt_ft - C.PLACES[tw[0]][2], 0) / 700.0, tw[1] / (gs * 1.852) * 60)
            return "lost low near place (landing estimated)", tw[0], tw[1], r.mdt + pd.Timedelta(minutes=m)
    return "lost en route", k, d, pd.NaT

pieces = []
for ic, g in T.sort_values("t").groupby("icao"):
    g = g.reset_index(drop=True)
    g["piece"] = ((g.t.diff() > 300) | (g.on_ground != g.on_ground.shift())).cumsum()
    parts = [q for _, q in g.groupby("piece")]
    for n, q in enumerate(parts):
        r0, r1 = q.iloc[0], q.iloc[-1]
        if q.on_ground.all():
            k, d = C.nearest_place(r0.lat, r0.lon)
            pieces.append(dict(icao=ic, reg=r0.reg, type=r0.type, kind="ground", start_mdt=r0.mdt, end_mdt=r1.mdt, n=len(q),
                               start_kind="on the ground", start_place=k, start_place_km=round(d, 1),
                               end_kind="on the ground", end_place=k, end_place_km=round(d, 1)))
            continue
        gb = n > 0 and parts[n - 1].on_ground.all() and (r0.t - parts[n - 1].t.iloc[-1]) < 600
        ga = n + 1 < len(parts) and parts[n + 1].on_ground.all() and (parts[n + 1].t.iloc[0] - r1.t) < 600
        s = classify_start(q, gb); e = classify_end(q, ga)
        kf = C.hav(q.lat.values, q.lon.values, *C.FIRE); i = int(np.argmin(kf))
        pieces.append(dict(icao=ic, reg=r0.reg, type=r0.type, kind="airborne", start_mdt=r0.mdt, end_mdt=r1.mdt, n=len(q),
                           callsigns=" ".join(sorted(set(q.callsign) - {""})),
                           start_alt_ft=r0.alt_ft, start_gs_kt=r0.gs_kt, start_track=r0.track_deg, start_lat=round(r0.lat, 4), start_lon=round(r0.lon, 4),
                           start_kind=s[0], start_place=s[1], start_place_km=round(s[2], 1), start_est_mdt=s[3],
                           end_alt_ft=r1.alt_ft, end_gs_kt=r1.gs_kt, end_track=r1.track_deg, end_lat=round(r1.lat, 4), end_lon=round(r1.lon, 4),
                           end_kind=e[0], end_place=e[1], end_place_km=round(e[2], 1), end_est_mdt=e[3],
                           min_km_fire=round(float(kf[i]), 1), min_km_fire_mdt=q.mdt.iloc[i], min_km_fire_alt_ft=q.alt_ft.iloc[i]))
P = pd.DataFrame(pieces).sort_values(["reg", "start_mdt"]).reset_index(drop=True)
# chain airborne pieces separated only by a coverage gap (next piece first heard en route, kept flying across the gap)
fid, flights = 0, []
P["flight"] = -1
for reg, g in P.groupby("reg", sort=False):
    prev = None
    for i in g.index:
        r = P.loc[i]
        if r.kind != "airborne":
            prev = None; continue
        join = False
        if prev is not None and r.start_kind == "first heard en route":
            pr = P.loc[prev]
            gap_min = (r.start_mdt - pr.end_mdt).total_seconds() / 60
            dist = float(C.hav(pr.end_lat, pr.end_lon, r.start_lat, r.start_lon))
            v = dist / max(gap_min, 0.1) * 60 / 1.852
            join = gap_min <= 45 and v >= 0.5 * np.nanmean([pr.end_gs_kt, r.start_gs_kt])
        if not join:
            fid += 1
        P.at[i, "flight"] = fid; prev = i
F = []
for f, g in P[P.kind == "airborne"].groupby("flight"):
    a, z = g.iloc[0], g.iloc[-1]
    q = T[(T.reg == a.reg) & T.mdt.between(a.start_mdt, z.end_mdt)]
    kf = C.hav(q.lat.values, q.lon.values, *C.FIRE); i = int(np.argmin(kf))
    F.append(dict(flight=f, reg=a.reg, type=a.type, callsigns=" ".join(sorted(set(" ".join(g.callsigns.fillna("")).split()))),
                  pieces=len(g), first_heard_mdt=a.start_mdt, last_heard_mdt=z.end_mdt,
                  origin_kind=a.start_kind, origin=a.start_place, origin_km=a.start_place_km, departure_mdt=a.start_est_mdt,
                  first_alt_ft=a.start_alt_ft, first_gs_kt=a.start_gs_kt,
                  destination_kind=z.end_kind, destination=z.end_place, destination_km=z.end_place_km, arrival_mdt=z.end_est_mdt,
                  last_alt_ft=z.end_alt_ft, last_gs_kt=z.end_gs_kt, last_lat=z.end_lat, last_lon=z.end_lon, last_track=z.end_track,
                  mean_gs_kt=round(float(q.gs_kt.mean()), 0), min_km_fire=round(float(kf[i]), 1), min_km_fire_mdt=q.mdt.iloc[i],
                  min_km_fire_alt_ft=q.alt_ft.iloc[i]))
FL = pd.DataFrame(F)
FL["operator_inferred_from_callsign"] = FL.callsigns.map(lambda s: C.operator_of(s.split())[0])

# ---------------------------------------------------------------- 4. Wainwright visits: flights ending / starting in the Wainwright area
WK = C.WAINWRIGHT_KEYS
vis = []
for _, f in FL.iterrows():
    if f.destination in WK and (f.destination_kind != "lost en route" or (f.destination_km <= 15 and f.last_alt_ft - 2230 <= 4500)):
        vis.append(dict(reg=f.reg, type=f.type, callsigns=f.callsigns, flight=f.flight, event="arrival",
                        last_or_first_heard_mdt=f.last_heard_mdt, alt_ft=f.last_alt_ft, km_to_place=f.destination_km,
                        nearest_wainwright_place=f.destination, basis=f.destination_kind, est_on_ground_mdt=f.arrival_mdt,
                        other_end=f"{f.origin_kind}: {f.origin} ({f.origin_km} km)",
                        other_end_mdt=f.departure_mdt if pd.notna(f.departure_mdt) else f.first_heard_mdt))
    if f.origin in WK and (f.origin_kind != "first heard en route" or (f.origin_km <= 15 and f.first_alt_ft - 2230 <= 4500)):
        vis.append(dict(reg=f.reg, type=f.type, callsigns=f.callsigns, flight=f.flight, event="departure",
                        last_or_first_heard_mdt=f.first_heard_mdt, alt_ft=f.first_alt_ft, km_to_place=f.origin_km,
                        nearest_wainwright_place=f.origin, basis=f.origin_kind, est_on_ground_mdt=f.departure_mdt,
                        other_end=f"{f.destination_kind}: {f.destination} ({f.destination_km} km)",
                        other_end_mdt=f.arrival_mdt if pd.notna(f.arrival_mdt) else f.last_heard_mdt))
V = pd.DataFrame(vis).sort_values("last_or_first_heard_mdt").reset_index(drop=True)

# ---------------------------------------------------------------- 5. aircraft summary (supplied file + archive context)
rows = []
from shapely.geometry import Point
import geopandas as gpd
burn = gpd.read_file(C.BURN).to_crs(3400).geometry.union_all()
from pyproj import Transformer
TR = Transformer.from_crs(4326, 3400, always_xy=True)
def km_to_burn(lat, lon):
    x, y = TR.transform(np.asarray(lon), np.asarray(lat))
    return np.array([burn.distance(Point(a, b)) / 1000 for a, b in zip(np.atleast_1d(x), np.atleast_1d(y))])
for h, g in S.sort_values("ts").groupby("hex"):
    g = g.copy()
    g["kf"] = C.hav(g.lat.values, g.lon.values, *C.FIRE)
    g["kb"] = km_to_burn(g.lat.values, g.lon.values)
    g["kw"] = C.hav(g.lat.values, g.lon.values, *C.PLACES["CYWV Wainwright Airport"][:2])
    g["kc"] = C.hav(g.lat.values, g.lon.values, *C.PLACES["CFB Wainwright airfield (CZWN stn)"][:2])
    cs = cs_by.get(h, [])
    op, code = C.operator_of(cs)
    r0, r1 = g.iloc[0], g.iloc[-1]
    ent = {"hex": h, "registration": r0.registration, "type": r0.type, "description": r0.description,
           "callsigns_in_archive": " ".join(cs), "operator_inferred_from_callsign": op,
           "positions_in_file": len(g), "first_mdt": r0.mdt, "first_lat": r0.lat, "first_lon": r0.lon, "first_alt_ft": r0.alt_ft,
           "last_mdt": r1.mdt, "last_lat": r1.lat, "last_lon": r1.lon, "last_alt_ft": r1.alt_ft,
           "alt_min_ft": int(g.alt_ft.min()), "alt_max_ft": int(g.alt_ft.max())}
    for tag, col in [("fire_centre", "kf"), ("burn_outline", "kb"), ("CYWV", "kw"), ("CFB_airfield_CZWN", "kc")]:
        i = g[col].idxmin()
        ent[f"closest_km_{tag}"] = round(float(g.at[i, col]), 2); ent[f"closest_mdt_{tag}"] = g.at[i, "mdt"]
        ent[f"closest_alt_ft_{tag}"] = int(g.at[i, "alt_ft"])
    ent["minutes_first_in_file_after_front_at_fire"] = round(C.mins_after_1500(r0.mdt) - C.FRONT_AT_FIRE, 1)
    rows.append(ent)
AS = pd.DataFrame(rows).sort_values("first_mdt").reset_index(drop=True)
role = {}
for r in AS.itertuples():
    pv = V[V.reg == r.registration]
    if not len(pv):
        fl = FL[(FL.reg == r.registration) & (FL.last_heard_mdt >= pd.Timestamp(r.first_mdt)) & (FL.first_heard_mdt <= pd.Timestamp(r.last_mdt))]
        role[r.registration] = "no arrival or departure in the Wainwright area; " + "; ".join(
            f"{f.origin} -> {f.destination} ({f.destination_kind})" for f in fl.itertuples())
    elif r.operator_inferred_from_callsign:
        role[r.registration] = "air ambulance operator by public call-sign prefix (inference); shuttles to/from Wainwright area"
    else:
        role[r.registration] = ("operator not identified from public call-sign data; flight pattern (arrivals/departures Wainwright area) "
                                "noted, purpose unknown" if len(pv) else "operator and purpose not identified")
AS["role_note"] = AS.registration.map(role)
AS["wainwright_arrivals"] = AS.registration.map(lambda r: int(((V.reg == r) & (V.event == "arrival")).sum()) if len(V) else 0)
AS["wainwright_departures"] = AS.registration.map(lambda r: int(((V.reg == r) & (V.event == "departure")).sum()) if len(V) else 0)
AS.to_csv(os.path.join(C.OUT, "aircraft_summary.csv"), index=False)

# ---------------------------------------------------------------- 6. response timeline
v5 = json.load(open(os.path.join(C.V5, "arrival_summary.json")))
def v5med(p): return C.mins_after_1500(pd.Timestamp("2026-10-17 " + v5["points"][p]["scar_limited_cell_weighted"]["clock"]["median"]))
REF = {"front_at_fire": C.FRONT_AT_FIRE, "v5_range2_median": v5med("Range 2 (RSO)"), "v5_range15_median": v5med("Range 15 (RSO)"),
       "v5_roadhead_median": v5med("Road head"), "dnd_just_before_1600": 60.0}
ev = []
def add(t, reg, what, detail=""):
    m = C.mins_after_1500(t)
    ev.append(dict(mdt=pd.Timestamp(t).strftime("%Y-%m-%d %H:%M:%S"), registration=reg,
                   callsigns=" ".join(cs_by.get(AS.loc[AS.registration == reg, "hex"].iloc[0], [])) if reg in AS.registration.values else "",
                   event=what, detail=detail,
                   min_after_front_at_fire=round(m - REF["front_at_fire"], 1), min_after_v5_range2_median=round(m - REF["v5_range2_median"], 1),
                   min_after_1600=round(m - 60.0, 1)))
add("2026-10-17 15:27:53", "", "MODEL: cold front at the fire (v5 front plane 15:27.9)")
for p, lab in [("Range 2 (RSO)", "Range 2"), ("Range 15 (RSO)", "Range 15"), ("Road head", "road head")]:
    c = v5["points"][p]["scar_limited_cell_weighted"]["clock"]
    add("2026-10-17 " + c["median"], "", f"MODEL: v5 weighted median fire arrival, {lab}", f"10-90%: {c['p10']}-{c['p90']}")
add("2026-10-17 16:00:00", "", "STATEMENT: DND, incident 'just before 4:00' local time (Oct. 22 briefing)")
g = T[(T.reg == "C-GSAB") & (T.mdt.between("2026-10-17 15:20", "2026-10-17 15:40"))]
if len(g):
    i = g.km_fire.idxmin()
    add(g.at[i, "mdt"], "C-GSAB", "closest approach to fire (transit, descending toward Lloydminster)",
        f"{g.at[i, 'km_fire']:.1f} km, {g.at[i, 'alt_baro']} ft")
for _, f in FL.iterrows():
    if f.last_heard_mdt < pd.Timestamp("2026-10-17 15:00"):
        continue
    if f.origin_kind != "first heard en route":
        add(f.departure_mdt, f.reg, f"take-off, {f.origin}" + (" (estimated)" if "estimated" in f.origin_kind else " (heard)"),
            f"first heard {f.first_heard_mdt:%H:%M:%S} at {f.first_alt_ft:.0f} ft, {f.origin_km} km from {f.origin}")
    else:
        add(f.first_heard_mdt, f.reg, f"first heard en route", f"{f.first_alt_ft:.0f} ft, {f.origin_km} km from {f.origin}")
    if f.destination_kind != "lost en route":
        add(f.arrival_mdt, f.reg, f"landing, {f.destination}" + (" (estimated)" if "estimated" in f.destination_kind else " (heard)"),
            f"last heard {f.last_heard_mdt:%H:%M:%S} at {f.last_alt_ft:.0f} ft, {f.destination_km} km from {f.destination}")
    else:
        add(f.last_heard_mdt, f.reg, "last heard en route", f"{f.last_alt_ft:.0f} ft, {f.destination_km} km from {f.destination}")
TL = pd.DataFrame(ev).sort_values("mdt").drop_duplicates().reset_index(drop=True)
TL.to_csv(os.path.join(C.OUT, "response_timeline.csv"), index=False)
P.to_csv(os.path.join(C.OUT, "flight_pieces.csv"), index=False)
FL.to_csv(os.path.join(C.OUT, "flights.csv"), index=False)
V.to_csv(os.path.join(C.OUT, "wainwright_visits.csv"), index=False)
print(FL[["reg", "callsigns", "first_heard_mdt", "origin_kind", "origin", "departure_mdt", "destination_kind", "destination", "arrival_mdt", "last_heard_mdt"]].to_string())

# ---------------------------------------------------------------- 7. earliest reaction -> timing bound -> likelihood over v5 runs
# Inbound response flights: flights whose destination is in the Wainwright area, heard after the front reached the fire.
FR = FL[FL.destination.isin(WK) & (FL.last_heard_mdt > pd.Timestamp("2026-10-17 15:28"))].copy()
FR["departure_or_first_heard_mdt"] = FR.departure_mdt.fillna(FR.first_heard_mdt)
FR["departure_basis"] = np.where(FR.origin_kind == "take-off (heard)", "heard from take-off",
                                 np.where(FR.origin_kind.str.contains("estimated"), "estimated from climb-out (departure field below coverage)",
                                          "first heard en route (departure earlier, unknown)"))
FR["origin_label"] = np.where(FR.origin_kind == "first heard en route", "first heard " + FR.origin_km.astype(str) + " km from " + FR.origin, FR.origin)
FR["transit_min"] = ((FR.arrival_mdt.fillna(FR.last_heard_mdt) - FR.departure_or_first_heard_mdt).dt.total_seconds() / 60).round(1)
FR = FR.sort_values("departure_or_first_heard_mdt").reset_index(drop=True)
FR.to_csv(os.path.join(C.OUT, "inbound_response_flights.csv"), index=False)
print(FR[["reg", "callsigns", "origin", "departure_basis", "departure_or_first_heard_mdt", "destination", "arrival_mdt", "transit_min", "mean_gs_kt"]].to_string())
star = FR[FR.callsigns.str.startswith("STAR") & (FR.departure_basis == "heard from take-off")].iloc[0]
earliest = FR.iloc[0]
t_launch_public = C.mins_after_1500(star.departure_or_first_heard_mdt)   # first inbound public air-ambulance flight, heard from take-off
t_launch_any = C.mins_after_1500(earliest.departure_or_first_heard_mdt)  # earliest inbound flight of any operator (may be estimated)
FR_first_on_ground = FR.arrival_mdt.min()

runs = pd.read_csv(os.path.join(C.V5, "ensemble_runs.csv"))
runs_sat = pd.read_csv(os.path.join(C.V5, "scar_satellogic", "ensemble_runs.csv"))[["run", "weight"]].rename(columns={"weight": "weight_satellogic"})
runs = runs.merge(runs_sat, on="run", how="left")
PTS = {"Range 2 (RSO)": "arr_scar|Range 2 (RSO)", "Range 15 (RSO)": "arr_scar|Range 15 (RSO)", "Road head": "arr_scar|Road head"}
from scipy.stats import norm
def air_like(t_arr, t_launch, min_lag, sigma=3.0):
    """P(first air-ambulance launch at t_launch | fire reaches the point at t_arr). With no public basis for the incident-to-launch lag
    beyond a minimum, the likelihood is flat for any t_arr <= t_launch - min_lag and falls off over ~sigma minutes after it."""
    return norm.cdf((t_launch - min_lag - np.asarray(t_arr, float)) / sigma)
def summ(x, w):
    ok = np.isfinite(x) & (w > 0)
    q = [C.C5.wquant(x[ok], w[ok], v) for v in (0.5, 0.1, 0.9)]
    sh = lambda a, b: float(w[ok & (x >= a) & (x < b)].sum() / w[ok].sum())
    return {"median": C.clock(q[0]), "p10": C.clock(q[1]), "p90": C.clock(q[2]),
            "share_3:30-3:45": round(sh(30, 45), 3), "share_3:45-4:00": round(sh(45, 60), 3), "share_after_4:00": round(sh(60, 1e9), 3)}
checks = {}
RW = runs[["run", "case", "front", "schedule", "fuel", "ros_mult", "spring_factor", "shape", "weighted", "weight", "weight_satellogic"] + list(PTS.values())].copy()
for lab, tl in [("public_operator_first_launch", t_launch_public), ("any_operator_first_launch", t_launch_any)]:
    for ml in (0, 10, 20, 30):
        # the bound applies to the incident, which is at or after the fire reaching the point the soldiers were at; using the
        # LATEST of the three points (most conservative) per run
        tmax = runs[list(PTS.values())].max(axis=1).values
        L = air_like(tmax, tl, ml)
        key = f"{lab}_minlag{ml}"
        res = {"launch_mdt": C.clock(tl), "min_lag_min": ml, "bound_mdt": C.clock(tl - ml),
               "v5_weight_retained_default": round(float((runs.weight * L).sum() / runs.weight.sum()), 4),
               "v5_weight_retained_satellogic": round(float((runs.weight_satellogic * L).sum() / runs.weight_satellogic.sum()), 4)}
        if ml == 10:
            RW[f"air_like_{lab}"] = L.round(5)
            for wcol in ("weight", "weight_satellogic"):
                w = (runs[wcol] * L).values; w = w / w.sum()
                RW[f"{wcol}_v51_{lab}"] = w
        checks[key] = res
RW.to_csv(os.path.join(C.OUT, "run_weights_air.csv"), index=False)
upd = {}
for p, col in PTS.items():
    x = runs[col].values.astype(float)
    upd[p] = {"v5_default": summ(x, runs.weight.values), "v5_satellogic": summ(x, runs.weight_satellogic.values),
              "v51_air_default": summ(x, RW["weight_v51_public_operator_first_launch"].values),
              "v51_air_satellogic": summ(x, RW["weight_satellogic_v51_public_operator_first_launch"].values)}
lags = {}
for lab, tm in [("front at fire 15:28", REF["front_at_fire"]), ("v5 Range 2 median", REF["v5_range2_median"]),
                ("v5 road head median", REF["v5_roadhead_median"]), ("v5 Range 15 p90 (latest 10-90% edge)",
                 C.mins_after_1500(pd.Timestamp("2026-10-17 " + v5["points"]["Range 15 (RSO)"]["scar_limited_cell_weighted"]["clock"]["p90"]))),
                ("DND 'just before 4:00' (taken as 15:55-16:00)", 57.5)]:
    lags[lab] = {"to_first_public_air_ambulance_take_off_min": round(t_launch_public - tm, 1),
                 "to_first_inbound_take_off_any_operator_min": round(t_launch_any - tm, 1),
                 "to_first_estimated_on_ground_in_Wainwright_area_min": round(C.mins_after_1500(FR_first_on_ground) - tm, 1)}
J["air_timing_check"] = {
    "first_inbound_public_air_ambulance": {k: (str(v) if isinstance(v, pd.Timestamp) else v) for k, v in star.items()},
    "earliest_inbound_any_operator": {k: (str(v) if isinstance(v, pd.Timestamp) else v) for k, v in earliest.items()},
    "logic": ["Incident (fire overtakes the position) -> casualties treated and moved by ground ('cared for first in Wainwright', "
              "DND Oct. 22) -> request for air transport -> crew launches -> transit -> landing in the Wainwright area.",
              "Only the launch and transit are seen. Every earlier step takes zero or more minutes, so the first inbound launch is an "
              "UPPER bound on the incident time: incident <= launch - (minimum lag).",
              "There is no public basis for the length of the ground phase, so the likelihood over incident times is flat below the "
              "bound. The air data therefore cannot favour any time inside the v5 window over another, or v5 over DND's 'just before 4:00'."],
    "bounds_and_v5_weight_retained": checks,
    "implied_lags": lags,
    "updated_arrival_summary": upd,
    "reading": "Every weighted v5 run reaches Range 2, Range 15 and the road head well before the bound, so v5.1's air-informed "
               "arrival windows are identical to v5's.",
}
print(json.dumps(J["air_timing_check"], indent=1, default=str)[:5000])

# ---------------------------------------------------------------- 8. burn / downwind proximity, approach paths, vertical irregularity
NEAR = T[T.km_fire <= 25].copy()
NEAR["km_burn"] = km_to_burn(NEAR.lat.values, NEAR.lon.values).round(2)
NEAR["bearing_from_fire"] = C.bearing(C.FIRE[0], C.FIRE[1], NEAR.lat.values, NEAR.lon.values).round(0)
NEAR["downwind_sector_100_150"] = NEAR.bearing_from_fire.between(100, 150)
NEAR["min_after_front_at_fire"] = ((NEAR.mdt - pd.Timestamp("2026-10-17 15:00")).dt.total_seconds() / 60 - C.FRONT_AT_FIRE).round(1)
px = NEAR[["reg", "type", "callsign", "mdt", "lat", "lon", "alt_baro", "gs_kt", "km_fire", "km_burn", "bearing_from_fire",
           "downwind_sector_100_150", "min_after_front_at_fire", "in_supplied_file"]].copy()
px["mdt"] = px.mdt.dt.strftime("%Y-%m-%d %H:%M:%S"); px["km_fire"] = px.km_fire.round(2)
px.to_csv(os.path.join(C.OUT, "proximity_burn_downwind.csv"), index=False)
ap = []
for _, v in V.iterrows():
    f = FL[FL.flight == v.flight].iloc[0]
    q = T[(T.reg == v.reg) & T.mdt.between(f.first_heard_mdt, f.last_heard_mdt) & (T.km_cywv <= 40)].copy()
    if len(q) < 3:
        continue
    # direct line between the far end of this 40 km portion and the Wainwright place nearest the near end
    far = q.iloc[0] if v.event == "arrival" else q.iloc[-1]
    wl = C.PLACES[v.nearest_wainwright_place]
    x, y = TR.transform(q.lon.values, q.lat.values); xf, yf = TR.transform(far.lon, far.lat); xw, yw = TR.transform(wl[1], wl[0])
    ux, uy = (xw - xf), (yw - yf); L = np.hypot(ux, uy)
    xt = ((x - xf) * uy - (y - yf) * ux) / L / 1000       # signed cross-track, km (positive = left of the direct line)
    kb = km_to_burn(q.lat.values, q.lon.values)
    brg = C.bearing(C.FIRE[0], C.FIRE[1], q.lat.values, q.lon.values)
    dw = (brg >= 100) & (brg <= 150) & (q.km_fire.values <= 25)
    ap.append(dict(reg=v.reg, callsigns=v.callsigns, event=v.event, heard_mdt=v.last_or_first_heard_mdt, n=len(q),
                   min_km_to_burn=round(float(kb.min()), 2), min_km_to_fire_centre=round(float(q.km_fire.min()), 2),
                   max_abs_cross_track_km=round(float(np.abs(xt).max()), 2), positions_in_downwind_sector=int(dw.sum()),
                   alt_range_ft=f"{q.alt_ft.min():.0f}-{q.alt_ft.max():.0f}",
                   direct_line_from=f"{far.lat:.3f},{far.lon:.3f}", direct_line_to=v.nearest_wainwright_place))
AP = pd.DataFrame(ap); AP.to_csv(os.path.join(C.OUT, "approach_paths.csv"), index=False)

def vert_stats(q):
    q = q.set_index("mdt").alt_ft.dropna()
    if len(q) < 6:
        return None
    r = q.resample("10s").mean().interpolate(limit=6).dropna()
    d = r.diff(3).dropna()                                  # 30 s altitude change
    sgn = np.sign(d[d.abs() >= 150])
    return dict(n=len(q), minutes=round((q.index[-1] - q.index[0]).total_seconds() / 60, 1),
                reversals_ge150ft=int((sgn.diff().abs() > 0).sum()),
                reversals_per_10min=round(float((sgn.diff().abs() > 0).sum()) / max((q.index[-1] - q.index[0]).total_seconds() / 600, 0.1), 2),
                alt_30s_change_sd_ft=round(float(d.std()), 0))
vs = []
for _, pc in P[P.kind == "airborne"].iterrows():
    q = T[(T.reg == pc.reg) & T.mdt.between(pc.start_mdt, pc.end_mdt)]
    for zone, keys in [("Wainwright area (<=40 km of CYWV)", WK), ("other airfield (<=40 km)", [k for k in C.PLACES if k not in WK])]:
        dmin = np.min([C.hav(q.lat.values, q.lon.values, C.PLACES[k][0], C.PLACES[k][1]) for k in keys], axis=0)
        z = q[(dmin <= 40) & (q.alt_ft <= 12000)]
        s = vert_stats(z)
        if s:
            vs.append(dict(reg=pc.reg, type=pc.type, piece_start=pc.start_mdt, zone=zone, **s,
                           before_2000=bool(pc.start_mdt < pd.Timestamp("2026-10-17 20:00"))))
VS = pd.DataFrame(vs); VS.to_csv(os.path.join(C.OUT, "vertical_irregularity.csv"), index=False)
vs_sum = VS.groupby("zone")[["reversals_per_10min", "alt_30s_change_sd_ft"]].median().round(2).to_dict()

# the C510 near the front
g = T[(T.reg == "C-GSAB") & T.mdt.between("2026-10-17 15:15", "2026-10-17 15:40")].copy()
g["front_min"] = C.front_time_at(g.lat.values, g.lon.values)
g["t_min"] = (g.mdt - pd.Timestamp("2026-10-17 15:00")).dt.total_seconds() / 60
g["min_after_surface_front"] = (g.t_min - g.front_min).round(1)
fp = g[["mdt", "lat", "lon", "alt_baro", "gs_kt", "track_deg", "vrate_fpm", "km_fire", "front_min", "min_after_surface_front"]].copy()
fp["mdt"] = fp.mdt.dt.strftime("%H:%M:%S"); fp["front_min"] = fp.front_min.round(1); fp["km_fire"] = fp.km_fire.round(1)
fp.to_csv(os.path.join(C.OUT, "front_passage_c510.csv"), index=False)
J["linkage"] = {
    "low_traffic_gap_within_80km_mdt": J["airspace_null_check"]["gap_below_12000ft_mdt"],
    "positions_within_25km": int(len(NEAR)), "first_position_within_25km_mdt": str(NEAR.mdt.min()),
    "positions_within_1km_of_burn": int((NEAR.km_burn <= 1).sum()),
    "min_km_to_burn_any_aircraft": float(NEAR.km_burn.min()) if len(NEAR) else None,
    "positions_in_downwind_sector": int(NEAR.downwind_sector_100_150.sum()),
    "downwind_sector_first_mdt": str(NEAR[NEAR.downwind_sector_100_150].mdt.min()) if NEAR.downwind_sector_100_150.any() else None,
    "approach_paths_max_cross_track_km_median": float(AP.max_abs_cross_track_km.median()) if len(AP) else None,
    "vertical_irregularity_median_by_zone": vs_sum,
    "c510_at_front": {"closest_km_fire": float(g.km_fire.min()), "alt_at_closest": str(g.loc[g.km_fire.idxmin(), "alt_baro"]),
                      "min_after_surface_front_range": [float(g.min_after_surface_front.min()), float(g.min_after_surface_front.max())]},
}
print(json.dumps(J["linkage"], indent=1, default=str))
json.dump(J["provenance"], open(os.path.join(C.OUT, "provenance_check.json"), "w"), indent=1, default=str)
json.dump({k: J[k] for k in ("air_timing_check", "airspace_null_check", "linkage")}, open(os.path.join(C.OUT, "air_timing_check.json"), "w"),
          indent=1, default=str)
print("done")
