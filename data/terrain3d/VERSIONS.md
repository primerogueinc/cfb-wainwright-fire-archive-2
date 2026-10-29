# 3D reconstruction videos: versions and the model behind each

The 3D videos follow the working-model versions. A model version bumps only when it is published, so every video made before Rubicon Ridge is published is a **working model v4** product, and lives in part 1 (frozen Oct. 27, 2026). Rubicon Ridge (**working model v5**) versions are built here, in part 2, and are **working results, not yet published** until v5 is.

## Model timeline (for the videos)

| When (2026, EDT) | Model | What it means for the videos |
|---|---|---|
| Oct. 24 | v3 published ("the noon photograph") | Start 3:30-7 a.m., 1.4-2.2 km east of Range 11. |
| Oct. 25 | v4 built (ensemble, smoke, scar fit, edge check) | Front at the fire 3:16 p.m. (3:25 allowance). Best-fit morning fire v4b, start ~5:00 a.m. Arrival: OSM "Range 2 point" 3:21, road head 3:26. Range positions from OpenStreetMap. |
| Oct. 26, 9 a.m. | 2024 Range Standing Orders added to v4 | Official grids: **Range 2 is ~375 m east of the OSM point** (188 m from the road head); **Range 15 (Field Firing) is the road head**; **Range 11 is ~240 m west of the OSM point**; **Range 13 sits ~930 m south of Range 11**, not where OSM puts it. v4 arrival at the official Range 2 grid 3:25, Range 15 3:26. |
| Oct. 26, 4:30 p.m. | Rubicon Ridge (v5) started | — |
| Oct. 27 | Part 1 frozen with v1-v4 as published | v5 work moves to part 2. |
| Oct. 27 | v5 working results | Front at the fire **~3:28 p.m.** (moving line fitted to 9 stations/sensors), not 3:16. Weighted median arrival: official Range 2 **3:39**, Range 15 **3:42**, road head **3:42**, v4's OSM point 3:36. Measured from the front, close to v4 (+6.3 vs +5.3 min at the OSM point). |

## Videos made so far (all working model v4; in part 1)

| Version | File (part 1, `data/terrain3d/`) | Commit | Made | Fire | Places | Notes |
|---|---|---|---|---|---|---|
| v4 3D-0 | `orbit_test/cfb-wainwright-3d-orbit-test.mp4` | f173cd2 | Oct. 26 00:35 | none (12:17 Landsat pixels, scar outline) | OSM | Height x3 (labelled). |
| v4 3D-0 | `pov/pov_sheet.jpg` and `pov/` | f173cd2 | Oct. 26 00:35 | 12:17 Landsat pixels | OSM (incl. OSM "Range 13") | Sightline counts. The "Range 13" view is at the OSM position, which the Standing Orders put elsewhere. |
| v4 3D-1 | `fire_draft/cfb-wainwright-3d-fire-draft1.mp4` | a9eb762 | Oct. 26 00:49 | **Simplified** spread (distance timing, v2/v3-style run 3:20-3:35), not v4's fire | OSM | Labelled v4 but its fire is not v4's. Superseded. Says "height x3", but true scale. |
| v4 3D-2 | `fire_draft/cfb-wainwright-3d-fire-draft2.mp4` | 50a8f4d | Oct. 26 02:51 | v4b morning fire from 8 a.m.; simplified run | OSM | Superseded. Says "height x3", but true scale. |
| v4 3D-3 | `fire_draft/cfb-wainwright-3d-timeline-draft3.mp4` | f425ffd | Oct. 26 03:47 | v4b morning + v4's best-fit post-front fire (315°, 8 km/h); front 3:16 | OSM | Full timeline 4:30 a.m.-3:45 p.m., escape scenarios. Says "height x3", but true scale. Superseded by 3D-5. |
| v4 3D-4 | `fire_draft/cfb-wainwright-3d-timeline-draft4.mp4` | 148d30e | Oct. 26 04:15 | as 3D-3 | OSM | Slowed key moments, lighter night. Same label error. Superseded by 3D-5. |
| **v4 3D-5** | `fire_draft/cfb-wainwright-3d-timeline-draft5.mp4` | d7f96f3 | Oct. 27 02:18 | as 3D-3 | OSM | **Current v4 timeline.** Scale label and figure placement corrected. |
| **v4 3D-E1** | `eye_level/cfb-wainwright-3d-eye-level-draft1.mp4` | d239e86 | Oct. 27 02:12 | as 3D-3; v4 smoke visibility | OSM Range 2 point, road head | **Current v4 eye-level.** "Range 2 point" is v4's OSM point; the official Range 2 grid is 375 m east, near the road head. |

SHA-256 (first 16 hex): 3D-0 orbit `6176f1fbd337de41`, 3D-1 `9ddaaec2e2dbbc3b`, 3D-2 `011afbabae28a501`, 3D-3 `7449a1cf754c8cea`, 3D-4 `047ad8018cdb88d2`, 3D-5 `fb1d9bcad8a12cc3`, E1 `c6762a32b9385660`.

## Videos made for v5 (Rubicon Ridge; working results, not yet published; in part 2)

| Version | File (part 2, `data/terrain3d/v5/`) | Made | Fire | Places | Notes |
|---|---|---|---|---|---|
| **v5 3D-T3** | `cfb-wainwright-3d-timeline-v5-draft3.mp4` | Oct. 28 | v5 weighted median; nothing drawn before the 12:17 detection; flames true scale (Byram) | 2024 Range Standing Orders | **Current timeline.** Opens at the 11:42 a.m. wind warning. Escape scenarios from the official Range 2 grid (+0/+5/+10 min after the front: 100% / 66% / 34% of v5's weight get clear). |
| **v5 3D-E3** | `cfb-wainwright-3d-eye-level-v5-draft3.mp4` | Oct. 28 | as T3; smoke timing carried over from v4 (v5 has no smoke model) | official Range 2 grid; road head (Range 15 entrance) | **Current eye level.** Fire arrives about 3:39 and 3:42. |
| **v5 3D-D1** | `cfb-wainwright-3d-drone-v5-draft1.mp4` | Oct. 28 | as T3 | as T3 | **New.** Hypothetical drone following the fire's edge, down to about 300 m over Range 2. |
| v5 3D-T1, E1 | `...-timeline-v5-draft1.mp4`, `...-eye-level-v5-draft1.mp4` | Oct. 28 | v5, but with v4's modelled morning fire from about 5 a.m. | as above | Superseded by T3 and E3. |

SHA-256 (first 16 hex): v5 3D-T3 `23cd4f0aa3fb29f9`, 3D-E3 `7acb08319b90f4c5`, 3D-D1 `03402778fc67c389`; v5 3D-T1 `68f810559ca97a35`, v5 3D-E1 `4c737195b249b34c`.

Details, method and the re-run check: `data/terrain3d/v5/ABOUT.md`. Code: `scripts/terrain3d_v5/`.

The v5 videos use v5's weighted-median fire, not a single "best" run: many runs tie on fit weight, and the median reproduces v5's published arrival medians.
