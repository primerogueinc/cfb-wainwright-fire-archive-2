# Rubicon Cloud, working model v5.1 (Oct. 27, 2026)

Model v5.1 adds aircraft positions (ADS-B) to Rubicon Ridge (working model v5, `data/rubicon-ridge/model_v5/`). It sits alongside v5 and **changes none of v5's numbers.** The scripts are in `scripts/rubicon_ridge/model_v5_1/`. This is a working model, not a final finding. All times are MDT (UTC−6) on Oct. 17, 2026, unless marked Oct. 18.

## Bottom line

- **The air data do not narrow when the fire overtook the soldiers.** They put a loose upper limit on it: **before about 16:55** (estimated) or **16:58** (heard). That is an hour or more after every time in play:
  - v5's arrival windows at Range 2, Range 15 and the road head (15:34–15:58);
  - DND's "just before 4:00";
  - the front reaching the fire at 15:28.
- **v5.1's arrival windows are therefore identical to v5's:** Range 2 15:39 (15:34–15:51), Range 15 15:42 (15:36–15:58), road head 15:42 (15:36–15:56). The Satellogic refit gives the same values.
- **No aircraft responded during the fire's intense phase, as far as ADS-B shows.**
  - After a transiting business jet left at 15:34, no aircraft at or below 12,000 ft was heard within 80 km of the fire until **17:20**, a gap of 106 minutes.
  - The first inbound air-ambulance launches were:
    - a King Air (**C-FSUG**) climbing out of the Lac La Biche area at about **16:55** (estimated);
    - a STARS helicopter (**C-FLCI, call sign STAR03**) lifting off at Edmonton International at **16:58** (heard).
  - Both reached the Wainwright area at about **17:35–17:37** (estimated).
- **That delay fits DND's account.** DND said casualties were "cared for first in Wainwright, then flown to hospital". Between the fire's arrival and the first air launch:
  - **about 79 min** at v5's Range 2 median;
  - **about 61 min** after "just before 4:00", taken as 15:55–16:00;
  - **90 min** after the front reached the fire.
  - The air data cannot say how that time divided between ground treatment, ground transport and the request for aircraft.
- **From 17:35 until about 04:00 the air response was a steady shuttle.**
  - In the Wainwright area it made 14 arrivals and 15 departures.
  - Fixed-wing aircraft flew from Wainwright Airport (CYWV) to Edmonton, Calgary and Saskatoon.
  - STARS helicopters flew to the University of Alberta Hospital site in Edmonton and the Foothills Medical Centre site in Calgary.
  - The operators are inferred from public call-sign prefixes:
    - STARS;
    - Can West Corporate Air Charters, an Alberta Health Services air-ambulance contractor;
    - the Saskatchewan Government Air Ambulance Service.
  - The operators of three aircraft are not identified.
- **No aircraft avoided anything that the data can show.**
  - No aircraft was heard near the fire while it was burning hard. GOES saw heat from 15:20 to 15:50.
  - In the evening, no aircraft came within 4.1 km of the burn outline.
  - Only one aircraft entered the downwind sector (bearing 100–150° from the fire): SLG2 at 18:19, 22 km from the burn, at 9,300–9,600 ft.
  - Approaches were not unusually bumpy. Altitude scatter near Wainwright was lower than near other airfields.
  - There is no basis for any plume-height or plume-avoidance inference.

## What the supplied file is

Prime Rogue supplied `raw/adsb_night_low/wainwright_fire_night_low_near_base_all_positions.csv` (SHA-256 `2fb3f95b…2602`).

- It is a **filtered subset of the adsb.lol open archive** (ODbL 1.0) that part 1 already holds.
- **Selection:** positions within 40 km of 52.80 N, 110.90 W, at or below 12,000 ft barometric altitude, from about 15:30 MDT Oct. 17 to 03:58 Oct. 18.
- **Match:** every row matches an archived position to within 1 s and 0.0006°.
- **"Base":** `km_from_base` measures distance from that round-number point. It is 4.0 km north of the fire centre and 2.8 km west of CYWV, not CZWN.

Full details are in `raw/adsb_night_low/README.md` and `provenance_check.json`. Because the file is a window, v5.1 reads the same aircraft's **whole flights** from part 1's original trace files (`context_tracks.csv`). That shows where each flight came from and went.

## Key times (MDT)

| Time | What | Basis |
|---|---|---|
| 15:23–15:34 | C-GSAB (Cessna Citation Mustang, C510) descends toward Lloydminster. It passes 29.5 km north of the fire at 15:29 at about 15,800 ft. | heard; transit, not a response |
| **15:27.9** | Cold front reaches the fire (v5 front plane) | model |
| 15:34–15:51 / 15:36–15:58 / 15:36–15:56 | v5 fire arrival, 10–90%: Range 2 / Range 15 / road head (medians 15:39 / 15:42 / 15:42) | model |
| "just before 4:00" | DND, Oct. 22 briefing | statement |
| **15:34–17:20** | **No aircraft at or below 12,000 ft heard within 80 km of the fire** (106 min). Airliners and other traffic passed at 17,900 ft and above. One of them, listed as "C-FHYY C172", flew at airline speed under call sign ACA339, which looks like a database mismatch. | heard (absence) |
| ~15:18 → **~16:55** | C-FSUG (King Air 200, BE20; operator not identified) descends toward the Lac La Biche area, then climbs out of it toward Wainwright | estimated (departure field below coverage) |
| **16:58:15** | **C-FLCI STAR03** (STARS H145 / EC145) lifts off at Edmonton International (first position 2,850 ft, 1.6 km from the field) | heard |
| 17:29 → ~17:35 | C-FSUG last heard 18.8 km from CYWV, descending through 6,275 ft; estimated landing | estimated |
| 17:32 → ~17:37 | STAR03 last heard 13.8 km west of CYWV, 4,600 ft, heading 112° toward CYWV and the town; estimated landing in the Wainwright area | estimated |
| 17:34 → ~18:24 | C-GSAH SLG2 (Saskatchewan Government Air Ambulance, inferred) from Saskatoon to CYWV | take-off heard; landing estimated |
| 17:36 → ~18:08 | C-GLUR CWA921 (Can West, inferred) from Edmonton to CYWV | take-off heard; landing estimated |
| ~17:57 → ~19:03 | C-GWRS (BE20; operator not identified) from the Fort McMurray area to CYWV | estimated |
| ~19:08 | First departure from Wainwright: SLG2 toward Saskatoon (lands 19:41) | estimated take-off; landing heard |
| ~19:43, ~20:49 | CWA921 and C-GWRS from CYWV to Edmonton | estimated take-off; landing heard |
| ~21:37–21:39 | STAR03 and STAR11 (C-FDRP, BK117) leave the Wainwright area. Both track ends are at the University of Alberta Hospital site in Edmonton (22:59 and 23:07). | estimated take-off; landing heard |
| to ~04:00 Oct. 18 | Shuttles continue: STAR01 from Calgary, CWA927, CWA935, SLG3 and others | `wainwright_visits.csv` |

`response_timeline.csv` has every event, with minutes after the front at the fire, after the v5 Range 2 median and after 16:00.

## Method

1. **Provenance** (`v51_air.py` §1). Each row is matched to part 1's adsb.lol trace files by ICAO address and time (±1 s). `km_from_base` is fitted to a point. The selection is rebuilt and checked row for row.
2. **Whole flights** (§2–3). Positions for the 11 aircraft are taken from noon Oct. 17 to 06:00 Oct. 18.
   - The track is split into continuously heard pieces wherever a gap exceeds 5 minutes. Pieces separated only by a coverage gap are joined into flights.
   - Each flight's ends are classified as follows:
     - **take-off or landing heard:** on the ground, or within about 2,000 ft of a nearby field's elevation;
     - **estimated:**
       - first or last heard low within 8 km of a field, heading suitably; or
       - climbing out of, or descending toward, a listed place along the ground track; or
       - lost low within 35 km of a place.
     - **en route:** neither of the above.
   - Estimated times assume 700 ft/min descents and 1,000–1,500 ft/min climbs, or the observed rate. They also allow for distance at the observed ground speed. Treat them as ±3–5 min.
   - Fixed-wing aircraft can only use aerodromes. Helicopters may also use the town or a hospital helipad. Place coordinates are approximate (`v51_common.py`).
3. **Operators** come only from public call-sign prefixes:
   - STARnn → STARS;
   - CWAnnn → Can West Corporate Air Charters (call sign CAN WEST);
   - SLGn → Saskatchewan Government Air Ambulance Service (LIFEGUARD).

   Sources are Nav Canada's list of Canadian ICAO designators and public airline lists. **These are inferences.** C-FSUG, C-GWRS and C-GSAB used their registrations as call signs; their operators and purposes are not identified here, and registry owners were not looked up. Nothing in v5.1 says who was aboard any aircraft.
4. **Air-informed timing check** (§7; `air_timing_check.json`, `run_weights_air.csv`).
   - The chain is: incident, then ground care and transport, then a request, then launch, then transit.
   - Only the last two steps are seen, so the first inbound launch is an **upper bound**: incident ≤ launch − minimum lag.
   - Each v5 run is weighted by L = Φ((T_launch − m − t)/3 min), where t is the latest of its Range 2, Range 15 and road-head arrivals.
   - The check is run for T_launch = 16:58 (STARS, heard) or 16:55 (any operator, estimated), and for m = 0, 10, 20 or 30 min.
   - In every combination, **100% of v5's weight (default and Satellogic) survives.** The v5.1 arrival summary, with m = 10 and STARS, equals v5's to the minute.
   - Below the bound the likelihood is flat. No public source gives the length of the ground phase, and assuming one would decide the answer.
5. **Avoidance and turbulence checks** (§8):
   - distance from each position within 25 km to the burn outline (Satellogic vectors) and the bearing from the fire (`proximity_burn_downwind.csv`);
   - cross-track deviation from a direct line on each Wainwright approach and departure (`approach_paths.csv`);
   - 30-second altitude changes and climb/descent reversals of 150 ft or more, near Wainwright and near other fields (`vertical_irregularity.csv`);
   - the C510's position relative to the v5 front plane (`front_passage_c510.csv`).
6. **Isochrones** (`v51_figures.py`). For each 40 m cell, the weighted median arrival over v5's 1,440 weighted scar-limited runs, contoured at 15:30–16:00 (`v5_isochrones_weighted_median.geojson`). This is model output.

## Model linkage, in detail

- **Front passage.** The only aircraft aloft nearby when the front reached the fire was C-GSAB. It was 29.5 km north of the fire at 15:29, at 15,500–16,500 ft, descending steadily at 1,100–1,700 ft/min with no altitude reversals.
  - At its position the surface front had passed about 20 minutes earlier, on the v5 plane.
  - Its profile shows nothing attributable to the front. It was in post-frontal air, well above the boundary layer and far from the fire. This is weak evidence either way.
- **Burn and downwind sector.** No aircraft position is within 4.1 km of the burn outline. The closest was SLG2 departing at 01:47 Oct. 18, 4.1 km north of the burn at 4,875 ft.
  - Every approach was flown in the evening or at night, long after the intense phase.
  - Large cross-track offsets on King Air arrivals (median 10 km) come from joining the circuit north and east of CYWV. They are not detours around the burn.
- **Turbulence.** No altitude reversals of 150 ft or more were found on any approach or departure within 40 km of CYWV (30 legs). The median 30-second altitude-change scatter was 84 ft near Wainwright, against 158 ft near other fields (55 legs). Nothing here suggests turbulence near the burn in the evening.
- **Plume height.** No inference is possible. There were no low aircraft in the downwind sector during the fire's active phase, and GOES saw no visible plume because of cloud.

## Limits

- **ADS-B coverage.** Volunteer receivers heard nothing near the base below about **4,250 ft barometric**, about 2,000 ft above the ground, so final approaches, landings and helicopter movements below that are not seen.
  - The supplied file also drops everything above 12,000 ft and everything more than 40 km from 52.80 N, 110.90 W. v5.1 reads whole flights from part 1 to get around this.
  - Short low hops, such as helicopter movements between the scene, the hospital and the airfields, can be entirely invisible.
- **Military aircraft are absent.** DND says RCAF helicopters carried casualties, but none appears. Military transponders often do not broadcast ADS-B, and only one aircraft flagged military appears anywhere in part 1's pull. **A missing track proves nothing.** An earlier response by unseen aircraft cannot be ruled out, and would tighten the bound.
- **Barometric altitude** is uncorrected for local pressure and is given in 25 ft steps. Times come from the receivers and are rounded to the second in the supplied file.
- **The bound is one-sided.** The air data can rule out an incident after about 16:55. They cannot favour 15:35 over 15:55. The lags quoted (61–90 min) are arithmetic, not findings about anyone's response.
- **Landing and take-off times** marked "estimated" are extrapolations of ±3–5 min. Places are inferred from where tracks start and end and are not confirmed.
- **Operators** are inferred from call-sign prefixes only.

## Files

| File | Contents |
|---|---|
| `provenance_check.json` | Supplied file vs the adsb.lol traces: row match, the "base" point, the reconstructed selection |
| `aircraft_summary.csv` | One row per aircraft in the file: first and last position; closest approach to the fire centre, the burn outline, CYWV and the CFB airfield (CZWN); call signs; inferred operator |
| `flights.csv`, `flight_pieces.csv` | Whole flights, noon Oct. 17 to 06:00 Oct. 18, with inferred origins, destinations and times |
| `wainwright_visits.csv` | Arrivals into and departures from the Wainwright area, with estimated landing and take-off times |
| `inbound_response_flights.csv` | Every flight into the Wainwright area after 15:28: departure basis, transit time, ground speed |
| `response_timeline.csv` | Chronological events with minutes after the front, the v5 Range 2 median and 16:00 |
| `airspace_1400_1800.csv` | All 36 aircraft heard within 80 km of the fire, 14:00–18:00 (the null check) |
| `air_timing_check.json` | The bound, the implied lags, v5 weight retained, the updated arrival summary and the linkage results |
| `run_weights_air.csv` | All 2,160 v5 runs: v5 weights, air likelihood and v5.1 weights (default and Satellogic) |
| `proximity_burn_downwind.csv`, `approach_paths.csv`, `vertical_irregularity.csv`, `front_passage_c510.csv` | Linkage checks |
| `context_tracks.csv` | Every position of the 11 aircraft, noon Oct. 17 to 06:00 Oct. 18, from part 1's adsb.lol traces (ODbL) |
| `v5_isochrones_weighted_median.geojson` | v5 weighted-median fire-arrival contours (model) |
| `fig_tracks_map.png` | Tracks coloured by time: whole flights, the supplied 40 km window, and a zoom with the burn and v5 isochrones |
| `fig_time_alt_dist.png` | Altitude and distance from the fire against time, per aircraft, with the front and v5 windows |
| `fig_response_timeline.png` | The air response as a Gantt chart, 3–6:45 p.m. and the whole night |
| `anim_tracks_1500-0400.gif` | 10-minute frames with 30-minute trails (0.3 MB) |

## Response-flight maps (`maps/`)

Plain maps of the same supplied positions (the 40 km, at-or-below-12,000-ft window), one colour per aircraft:

| File | What |
|---|---|
| `maps/aircraft_tracks_map.png` | All 11 aircraft, full tracks in the window, with airfields and heliports (Wainwright Airport CYWV runway 13/31; Field 21 / CZWN; Camp Wainwright helipad CFF7; Wainwright Health Centre heliport CWH2), each with its distance from the fire centre, and 5–40 km distance rings |
| `maps/aircraft_tracks_animation.mp4` | The same as an animation, 2-minute steps, 3:31 p.m. Oct. 17 to 3:58 a.m. Oct. 18 MDT, 20-minute trails, altitude labels (31 s) |
| `maps/aircraft_range_topo.png` | Tracks over the ranges and terrain (NRCan MRDEM-30 hillshade, 10 m contours), with official range grids (2024 Range Standing Orders), OpenStreetMap range areas and danger-area boundaries, 2–10 km rings, and the time each aircraft dropped below or climbed back into receiver coverage |

Distances are straight-line from the fire centre (52.7645 N, 110.9052 W): Wainwright Airport about 4.7 km, Camp Wainwright helipad about 7 km, Health Centre heliport about 9 km, Field 21 about 14.7 km. Airfield positions: Transport Canada data via SkyVector and AirportRoutes. Receivers heard nothing near the base below about 4,250 ft (about 2,000 ft above ground), so the maps cannot show landings. Code: `scripts/rubicon_ridge/model_v5_1/maps/` (`python3 air_anim.py <positions.csv> <out_dir>`, `python3 air_range_topo.py <positions.csv> <out.png>`).

## Reproduce

```
cd scripts/rubicon_ridge/model_v5_1
python3 v51_air.py            # needs part 1 next to this repository, or CFB_ARCHIVE_1=/path
python3 v51_figures.py --gif  # isochrones need the stored v5 grids (V5_SCRATCH/grids); otherwise the committed GeoJSON is reused
```

## Credits and licences

- Aircraft positions: [adsb.lol](https://adsb.lol) open archive, © adsb.lol contributors, [ODbL 1.0](https://opendatacommons.org/licenses/odbl/1-0/). The derived tables here (`context_tracks.csv` and the others) are adapted from that database and carry the same licence.
- Burn outline: Prime Rogue's vector measurements from a Satellogic 70 cm image via SkyFi, Oct. 20, 2026. No imagery is stored or shown.
- Model, analysis and figures: Prime Rogue Inc., CC BY-NC 4.0.
