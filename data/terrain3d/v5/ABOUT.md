# 3D reconstruction videos, Rubicon Ridge (working model v5)

**Working results, not yet published.** Rubicon Ridge is working model v5 (`data/rubicon-ridge/model_v5/`). Until it is published, these videos are drafts. The v4 versions stay in part 1 (`data/terrain3d/`). Prime Rogue Inc., CC BY-NC 4.0. Made Oct. 28, 2026.

Code: `scripts/terrain3d_v5/`, copied from part 1's `scripts/terrain3d/` and changed as below. The terrain, imagery and camera renders are the v4 ones (`scene.py`, true scale).

## Videos

| Version | File | What it shows |
|---|---|---|
| v5 3D-T1 | `cfb-wainwright-3d-timeline-v5-draft1.mp4` | The full day from the 4:30 a.m. forecast to the fire's run southeast, from above, with the escape scenarios |
| v5 3D-E1 | `cfb-wainwright-3d-eye-level-v5-draft1.mp4` | Eye level (1.7 m) from the official Range 2 grid (top) and the road head / Range 15 entrance (bottom), 11:30 a.m. until the modelled fire arrives |

Both videos are in this folder.

## What changed from the v4 videos

| | v4 videos (part 1) | v5 videos |
|---|---|---|
| Front at the fire | 3:16 p.m. (the CZWN gust) | **3:28 p.m.** (15:27.9; moving line fitted to 9 stations and sensors, toward 114.5°, 77 km/h). The 3:16 gust is shown at CZWN, 14.7 km west-northwest, and the fitted front is drawn crossing the scene |
| Fire | v4 case v4b and its best post-front run (one run) | **Weighted median of v5's 1,440 weighted runs**, cell by cell (`layers/v5_arrival_weighted_40m.tif`: median, 10th and 90th percentiles, burn probability). Before the front this is the median of v4's four morning cases (start 4:00-6:30 a.m.) as carried into v5 |
| Places | OpenStreetMap Range 2 point, Range 11 | **2024 Range Standing Orders grids:** Range 2, Range 11; the road head labelled "Range 15 / road head" (the Range 15 barrier is 16 m away) |
| Fire reaches Range 2 / road head | 3:21 / 3:25 (Range 2 point, OSM) | **3:39 / 3:42** (official Range 2 grid / road head; v5 weighted medians, 10-90%: 3:34-3:51 and 3:36-3:56) |
| Figures | moving until 3:16 | moving until the front reaches the fire (3:28) |
| Smoke at eye level | v4 smoke model | v5 has no smoke model. v4's road-head timing is carried over relative to the front and the fire: clear until the front, under 1 km 2 minutes after it, under 100 m 2 minutes before the fire (log-interpolated) |

## Re-running v5

The v5 arrival grids are not in the repository. `v5_rerun_grids.py` re-runs the 1,440 weighted scar-limited solves exactly as `v5_ensemble.py` (published confinement: revised scar + 150 m). Every run's arrival at the official Range 2 grid, Range 15, the road head and v4's Range 2 point matches the published `ensemble_runs.csv` to within 0.005 minutes. `v5_median_grid.py` builds the weighted per-cell quantiles. The weighted median reproduces the published medians: Range 2 15:39, Range 15 15:42, road head 15:42, v4's point 15:36.

A per-cell median is not itself one fire. It is the most likely arrival time at each place across the runs, and the video draws it as one fire. Its burned area (13.3 km² at 10 m) includes the 150 m margin the scar-limited runs may burn. About 1.0 km² of the mapped scar is not reached in at least half the weight; this is shown grey at the end as "timing unknown".

## Escape scenarios (`layers/scenarios_v5.json`, `scenarios_v5.py`)

The method follows v4 part 4 (part 1, `data/authentic-stoic-stress-test/part4_escape/`), simplified:

- **Routes and pace:** straight routes in 72 directions from the official Range 2 grid, jogging in kit at 8 km/h on flat ground.
- **Safe ground:** ground the run never reaches (outside the scar + 150 m), or ground that burned at least 5 minutes earlier.
- **Margin:** the fire must stay at least 3 minutes behind at every point.
- **Departure times:** measured from each run's own front time at the fire.

| Leaving after the front reaches the fire | Range 2 (official grid) | Road head |
|---|---|---|
| +0 min (about 3:28) | 100% | 100% |
| +5 min (about 3:33) | 66% | 87% |
| +10 min (about 3:38) | 34% | 46% |

Shares are of the v5 fit weight. On the median fire, the best route from Range 2 runs east (080°) and reaches ground the fire never burns after about 200 m. Official Range 2 lies 10 m outside the mapped scar edge.

These numbers are **more favourable than v4's** (80% / 43% / 10% from v4's Range 2 point), because the official Range 2 grid sits at the edge of the burn while v4's point was inside it. They share v4's limits: no smoke, embers or spotting, perfect knowledge of the direction, no time to react, and a smooth elliptical fire. Treat them as best-case upper bounds.

## Files

- `layers/points.json`: scene points. `pov` and `extra` hold the Standing Orders grids. `camera_ref` holds the v4 OpenStreetMap points the camera shots were framed on.
- `layers/modelled_arrival_minutes_v5.tif`: the 10 m arrival field used in the videos, in minutes after midnight MDT.
- `layers/v5_arrival_weighted_40m.tif`: the v5 model grid (40 m, EPSG:3400) with four bands: weighted median (minutes after 3:00 p.m.), p10, p90, and burn probability.
- `layers/scenarios_v5.json`, `layers/fire_meta.json`, `layers/povcam_*.json`: the escape scenarios, the fire metadata and the eye-level cameras.

## Illustrative

These parts are illustrative, not modelled:

- the AGL firing and every figure (representative, not to scale; DND has not said where C Company was);
- flame heights and the smoke column;
- the smoke timing at eye level.
