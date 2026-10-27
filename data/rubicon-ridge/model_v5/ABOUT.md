# Rubicon Ridge, working model v5 (Oct. 27, 2026)

This model reconstructs how the Oct. 17, 2026 grass fire at CFB Wainwright spread after the cold front arrived. It follows working model v4 (`data/authentic-stoic-stress-test/model_v4/`), which stays as published. Scripts are in `scripts/rubicon_ridge/model_v5/`. It is a working model, not a final finding.

## In plain English

- **The front reached the fire later than v4 assumed.** Nine weather stations and air sensors put the cold front over the fire at about **3:28 p.m.** (between 3:27 and 3:29 in leave-one-out checks). v4 used 3:16 p.m. (the airfield station, about 15 km west) and allowed for 3:25.
- **After the front, the fire reached the official Range 2 grid at about 3:39 p.m.** (weighted median; 10-90% range 3:34-3:51). It reached **Range 15 at about 3:42** (3:36-3:58) and **the road head at about 3:42** (3:36-3:56). On the model, no point burns before 3:30.
- **Measured from the front's arrival, the timing is close to v4's.** At v4's Range 2 point the fire arrives +6.3 min after the front (v4: +5.3) and at the road head +12.0 min (v4: +10.1). Most of the clock-time shift comes from the later front.
- **How fast the fire ran after the front is the largest unknown.** The final burn scar cannot settle it, because a fast fire and a slow fire leave the same shape. At the speed an independent expert said to expect (about 12.7 km/h), the fire reaches Range 2 at about 3:35 and the road head at about 3:37. At half the standard rate (about 3.5-5 km/h) it gets there at about 3:48-3:53.
- **The April 2026 planned burn did not clearly slow the fire.** Giving that ground a lower spread rate makes the fit to the scar slightly *worse*, not better.

## What changed from v4

| | v4 (Oct. 25) | v5 Rubicon Ridge |
|---|---|---|
| Front at the fire | 3:16 (CZWN) or 3:25 (allowance), same moment everywhere | Moving line fitted to 9 station/sensor times: toward 114.5°, 77 km/h, at the fire 15:27.9 (leave-one-out 15:27.5-15:28.5). Runs use 15:27 / 15:29 / 15:31 |
| Pre-front weather after 15:00 | Interpolated toward the 16:00 record (post-front air leaked in) | Held at the 15:00 observation until the front reaches each cell |
| Post-front wind | 5 schedules, including the slow "authentic stoic" veer | 5 schedules from **sustained** 10 m wind (never gusts): CZWN observed (61 km/h first half hour from the hourly vector mean, then 48 → 45 → 37 km/h, veering 290→307°); fixed 300°, 310° (Hughenden), 315° (v4's best) at 48 km/h then the CZWN hourly decline; RAP/HRDPS model mean at the fire |
| Post-front spread rate | Fixed 4.4 / 8 / 12 km/h | Canadian FBP grass, O-1a (matted) or O-1b (standing), 95% curing, ISI from hourly FFMC and sustained wind, × 0.5 / 1.0 / 1.25: 3.5-12.9 km/h head spread over the first 30 min |
| April 20-21, 2026 planned burn | Ordinary grass | Separate fuel class at 0.5× or 0.75× spread (1.0× run as an unweighted control) |
| Reference points | OSM "Range 2 point", road head | Official Range 2 and Range 15 grids (2024 Range Standing Orders), road head, plus v4's point for comparison |
| Runs | 240 | **2,160** (full balanced design, no reduction needed) |
| Best overlap with scar (Jaccard) | 0.51-0.55 | 0.536 (weighted runs); 0.548 including the control |

`fig_v4_vs_v5.png` and `v4_vs_v5.csv` give the side-by-side table.

## Method

1. **Front** (`v5_inputs.py`, `front_plane.json`): least-squares plane `t = a + b·x + c·y` through the front's arrival time at CZWN, four ACIS stations (Albert Hall, Bellshill, Gilt Edge North, Hughenden), two PurpleAir sensors (Wainwright, Castor), Coronation and Lloydminster. RMS residual 3.0 min; Bellshill (−7.1) and Castor (+4.4) are the outliers. Each cell meets the front at its own time; the post-front wind pattern travels with the front.
2. **Pre-front fire** (`v5_progression.py`): v4's four bracketing cases (start 4:00-6:30 a.m., morning spread 1.5-3×, width 1.5-3×) run to 15:35 with the 15:00 weather held. The outline of the 12:17 Landsat and 2:49 p.m. VIIRS detections is always counted as burned.
3. **Post-front spread** (`v5_common.py: spread_v5`): v4's 40 m Huygens/Dijkstra spread with an elliptical fire (ignition at the rear focus), now time-dependent. Shape is "steady" (FBP grass L:B = 1.1·WS^0.464, about 6.4 at 48 km/h) or "gusty" (L:B 3).
4. **Spring burn factors.** About six months of regrowth on burned ground, with no thatch or litter, carries less fuel. FBP O-1 has no load term, so a spread factor stands in for it. **0.75** is roughly the reduction grassland models apply for grazed or cut pasture, and **0.5** is roughly that for eaten-out pasture (CSIRO grassland model, Cheney, Gould & Catchpole 1998; values approximate). The weighting uses only these two values, as designed; 1.0 is a control.
5. **Ensemble** (`v5_ensemble.py`): 4 cases × 3 front times × 5 winds × 2 fuels × 3 multipliers × 3 spring factors × 2 shapes = 2,160 runs. Each run is solved twice:
   - **free**: anything but open water burns. The run is stopped when its area equals the scar's, and its overlap with the scar is the fit.
   - **scar-limited**: confined to the scar plus 150 m, as v4. This version supplies the arrival times and the GOES check.
6. **Fit** (`v5_fit.py`): as v4, w = exp((J − J_max)/0.05) over the 1,440 runs with spring factor 0.5 or 0.75. Two sensitivities:
   - the control runs weighted too;
   - weights multiplied by a GOES-18 term exp(−½(ln(model/obs)/ln 3)²) over the 15:10, 15:30 and 15:40 scans.
7. **Arrival times** at the cell containing each point, and at the earliest cell within 60 m. Official Range 2 lies 10 m and Range 15 18 m outside the mapped edge, which puts both on the edge at 40 m resolution. Free runs give the same medians.

## Results (`arrival_summary.json`, `ensemble_runs.csv`, `fit_weights.csv`)

**Arrival (weighted by fit to the scar, scar-limited runs; clock MDT):**

| Point | Median | 10-90% | Minutes after front at the fire | 3:30-3:45 / 3:45-4:00 / after 4:00 |
|---|---|---|---|---|
| Range 2 (official grid) | 15:39 | 15:34-15:51 | +9.7 (6.0-20.8) | 74% / 25% / 1% |
| Range 15 (official grid) | 15:42 | 15:36-15:58 | +13.3 (8.1-29.8) | 62% / 31% / 7% |
| Road head | 15:42 | 15:36-15:56 | +12.0 (7.8-27.3) | 65% / 30% / 5% |
| Range 2 point (OSM, v4) | 15:36 | 15:32-15:43 | +6.3 (3.8-14.5) | 93% / 6% / 0% |

Within 60 m, each point is reached about a minute earlier. Each 2-minute step in the front time moves every arrival by 2 minutes; minutes-after-front are unchanged. Every weighted run burns all three points.

**By spread rate (weighted medians):**

| Multiplier | Head spread (first 30 min) | Range 2 | Range 15 | Road head |
|---|---|---|---|---|
| 0.5 | 3.5-5.1 km/h | 15:48 (+19.1) | 15:54 (+24.1) | 15:53 (+23.5) |
| 1.0 | 7.0-10.3 km/h | 15:38 (+9.2) | 15:41 (+11.6) | 15:40 (+11.3) |
| 1.25 | 8.8-12.9 km/h | 15:36 (+7.2) | 15:39 (+9.1) | 15:38 (+8.8) |
| runs ≥ 11 km/h (near the expected ~12.7) | 11-12.9 km/h | 15:35 (15:32-15:38) | 15:37 (15:34-15:40) | 15:37 (15:33-15:40) |

**Best-fit parameters** (`fig_best_runs_vs_scar.png`). The best run is case v4c, front 15:27, fixed 315°, O-1b × 1.25, spring factor 0.75, steady shape (J = 0.536). Weight shares by factor:

- **Wind:** fixed 315° 30%, 310° 25%, CZWN observed veer 15%, 300° 15%, RAP/HRDPS 15%. The scar again favours north-westerly winds, as in v4.
- **Shape:** steady 58%, gusty 42%.
- **Spread rate, fuel type, pre-front case and front time: each about evenly split.** The scar does not constrain them.

**Implied post-front spread rate.** The weighted median is 7.7 km/h (10-90%: 3.5-11.6). This is effectively the prior range: the area-matched fit cannot tell speeds apart. For reference, FBP O-1a at most about 11.4 km/h, O-1b about 15 km/h, and an independent expert's expected value is about 12.7 km/h. Only 16% of the weight is on runs of 11 km/h or more.

**April 20 planned burn:**

- On the scar, 48.6% of the April burn burned again. The weighted runs burn 59% of it (10-90%: 54-62%).
- Lowering its spread rate *lowers* the best overlap: 0.548 at 1.0, 0.536 at 0.75, 0.519 at 0.5.
- Weighting the control runs too (sensitivity) gives 36% / 34% / 31% to factors 1.0 / 0.75 / 0.5. Arrivals then move 1-2 min earlier: Range 2 15:38, Range 15 15:41, road head 15:40.
- **Reading:** the scar does not support a large regrowth effect at 40 m. The model's misfit lies elsewhere (below).

**GOES-18 check** (indicative only; `goes_check` in the JSON):

- **No large growth before the front.** Modelled fire power at 15:10 is 0 MW for cases v4a, v4b and v4d and about 100 MW for v4c, against 76.5 MW observed. Pre-front growth from 15:10 to the front is 0.002-0.13 km².
- **After the front the model is hotter than GOES.** Medians are about 3,600 MW at 15:30 (observed 1,215) and about 2,000 MW at 15:40 (observed 336).
- **This check is crude.** In the model the whole perimeter of the roughly 2 km² fire ignites at once when the front passes. Fuel load and radiative fraction are assumed, and GOES read 0 MW at 15:20.
- **GOES-weighted sensitivity:** it favours the slowest runs (multiplier 0.5: 57%) and case v4c (99.9%, the only case with 15:10 fire power). Arrivals become Range 2 15:44 (+14.1 min), Range 15 15:48, road head 15:47. This bounds the late side; the overlap weighting stays the headline.

## Limitations

- **The fit cannot pin down speed.** Free runs are stopped at the scar's area, so the overlap measures shape only. Spread rate drives the arrival spread (road head +8.8 to +23.5 min across multipliers), and nothing in v5 constrains it except the GOES check, which is crude.
- **Overall fit is modest (J ≤ 0.55).** The model burns too much to the north-east, into the April burn and the 33-35 ha unburned pocket by the road head (every run burns the pocket, as in v4). It burns too little in the scar's south-east lobe beyond the boundary road and at its western tip. Something local is missing, such as tracks, fuel breaks or suppression.
- **Front time** is a plane fitted to 9 points with a 3-minute RMS residual; there is no ground observation at the fire.
- **Wind:**
  - The CZWN first-segment speed (61 km/h) is a derived sustained mean (hourly vector mean minus the pre-front minutes), higher than the 45-50 km/h ten-minute means.
  - All schedules are point winds applied uniformly across the grid.
- **Fuel:**
  - Curing is fixed at 95%, and FFMC comes from the airfield's hourly weather.
  - The spring-burn factors are judgement values (see method).
- **Resolution:**
  - The grid is 40 m, and the points lie 10-20 m from the mapped edge.
  - The scar is the Sentinel-2-derived `scar_patches_revised.geojson` (largest patch).
- **Satellogic outline:** the sharper 70 cm outline (`data/rubicon-ridge/satellogic_2026-10-20/`) did not exist when this was run. Refit command below.
- **Engine bug, found and fixed before any reported run.** The first v5 attempt fed float32 pre-front times into the spread engine. Rounding then let cells be "improved" endlessly: memory blew up and arrival times came out wrong. Its one finished block was discarded. `spread_v5` now forces float64 and uses an array heap, and its output was checked against a constant-wind case that reproduces v4's 4.7-5.2 min to v4's Range 2 point.

## Files

- `config.json`, `front_plane.json`, `wind_schedules.json`, `schedules.csv`: inputs.
- `prefront/`, `prefront_progression.csv`: the four pre-front cases to 15:35.
- `ensemble_runs.csv`: all 2,160 runs, with factors, fit, arrivals, GOES numbers and weights.
- `fit_weights.csv`: runs ranked by weight.
- `arrival_summary.json`: weighted and sensitivity arrival distributions at Range 2, Range 11, Range 15, the road head and v4's Range 2 point; weight shares; implied spread rate; GOES check.
- `fig_best_runs_vs_scar.png`, `fig_arrival_distributions.png`, `fig_v4_vs_v5.png`, `v4_vs_v5.csv`.

## Reproduce

```
cd scripts/rubicon_ridge/model_v5
python3 -I v5_inputs.py                 # front plane, wind schedules
python3 v5_progression.py ../../..      # pre-front cases
python3 v5_ensemble.py [--cases A,B]    # 2,160 runs; checkpoints per (case, front) block in $V5_SCRATCH/grids, restart skips done blocks
python3 v5_fit.py                       # weights, summaries, figures
# refit to another outline from the stored grids (about 15 s, no re-run):
python3 v5_fit.py --scar ../../../data/rubicon-ridge/satellogic_2026-10-20/<outline>.geojson --scar-mode union --tag satellogic
```

The arrival grids (float16, compressed, about 110 MB) live in the scratch folder and are not part of the repository. `V5_SCRATCH` sets their location. A refit with `--scar` keeps the scar-limited arrival runs on the original scar + 150 m; re-run `v5_ensemble.py --scar ... --tag X` to rebuild those too (about 10 min on 2 cores), then `v5_fit.py --tag X --grids $V5_SCRATCH/grids_X`.
