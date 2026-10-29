# CFB Wainwright Fire Archive, part 2 (2026)

Part 2 of Prime Rogue Inc.'s open archive of data and analysis on the deadly fire at Canadian Forces Base Wainwright, Alberta, on **October 17, 2026**. Maintained by [Prime Rogue Inc.](https://primerogueinc.com) alongside its [live coverage](https://primerogueinc.com/blog/cfb-wainwright-fire-live-updates/).

**Part 1** is [primerogueinc/cfb-wainwright-fire-archive](https://github.com/primerogueinc/cfb-wainwright-fire-archive). It holds the first evidence packages, working models v1–v4 as published (tag [`v4-published`](https://github.com/primerogueinc/cfb-wainwright-fire-archive/releases/tag/v4-published)), the earlier fires, ATIP records and the automated FIRMS and Sentinel-2 collection. Part 1 is now frozen. New work goes here.

Part 2 covers **Rubicon Ridge**, working model v5. It started on Oct. 26, 2026. The `data/rubicon-ridge/` and `scripts/rubicon_ridge/` folders were copied from part 1 as they stood on Oct. 27, so their history before that date is in part 1.

## What's here

| Path | Contents |
|---|---|
| `data/rubicon-ridge/DATA_SOURCES.md` | Rubicon Ridge's inputs: what each is for and where it is archived. |
| `data/rubicon-ridge/raw/` | Raw pulls for Oct. 16–18, kept unaltered. FIRMS (all sensors), PurpleAir 2-minute histories and size fractions, GOES-18/19 fire extracts, RAP, ECCC HRDPS, ECCC radar, 38 weather stations (SWOB), Alberta air quality, Alberta ACIS hourly stations, and Sentinel-3 SLSTR FRP. `README_pull.md` has the front-passage table. No API keys are stored. |
| `data/rubicon-ridge/derived/` | Hourly maximum-gust map and table (38 stations, 12–5 p.m. MDT). PurpleAir PM1/PM2.5/PM10 by hourly tranche: the 3:28 p.m. spike is coarse dust, not smoke. The front-arrival map (fire ~3:29 p.m. MDT). |
| `data/rubicon-ridge/derived/goes_visible/` | GOES-18/19 visible imagery, 2:30–5:30 p.m. MDT. **No identifiable smoke plume** (frontal cloud). 3.9 µm heat signal 3:20–3:50 p.m. |
| `data/rubicon-ridge/derived/polar_visible/` | Every VIIRS and MODIS pass, 1–6 p.m. MDT, with geometry and cloud. No polar pass imaged the intense phase. The leads are Aqua MODIS at 3:18 p.m. and NOAA-20 at 2:50 p.m. The imagery itself has not been pulled yet. |
| `data/rubicon-ridge/model_v5/` | Working model v5: 2,160-run ensemble, scar-fit weights, arrival times at Range 2, Range 15 and the road head, and the v4-vs-v5 comparison. `scar_satellogic/` refits it to the 70 cm burn outline. **Working results, not yet published.** See `ABOUT.md`. |
| `data/rubicon-ridge/raw/adsb_night_low/` | Supplied ADS-B positions near the base, 3:31 p.m. Oct. 17 to 3:58 a.m. Oct. 18, unaltered, with checksum. A filtered subset of the adsb.lol archive. See its `README.md`. |
| `data/rubicon-ridge/model_v5_1_rubicon_cloud/` | Working model v5.1 "Rubicon Cloud": the air-response timeline from ADS-B and how it links to v5. No aircraft at or below 12,000 ft was heard within 80 km of the fire from 3:34 to 5:20 p.m. The first inbound air-ambulance launches were at about 4:55 (estimated) and 4:58 p.m. (STARS, heard). This is an upper bound only, so v5's arrival windows are unchanged. **Working results, not yet published.** See `ABOUT.md`. |
| `data/rubicon-ridge/satellogic_2026-10-20/` | Prime Rogue's measurements from a Satellogic 70 cm image taken Oct. 20, 2026, 12:56 p.m. MDT. Includes the burn outline, unburned islands, edge runs along tracks, ground classes, objects near the ranges, distances from Range 2, Range 15 and the road head to the burn edge, and the image footprint. **Vectors and tables only. No imagery** (see Licences). See `ABOUT.md`. |
| `data/terrain3d/VERSIONS.md` | Which 3D reconstruction video goes with which working model: the v4 videos in part 1 (frozen) and the Rubicon Ridge (v5) versions. |
| `data/terrain3d/v5/` | Rubicon Ridge (v5) 3D videos, working results: full timeline, eye level from the official Range 2 grid and the road head, and a hypothetical drone following the fire. Code: `scripts/terrain3d_v5/`. |
| `scripts/rubicon_ridge/` | The code behind all of the above. |

## Running model v5

Model v5 reuses the v4 grid, inputs and code from part 1. Check out both repositories side by side:

```
git clone https://github.com/primerogueinc/cfb-wainwright-fire-archive
git clone https://github.com/primerogueinc/cfb-wainwright-fire-archive-2
```

If part 1 is somewhere else, set `CFB_ARCHIVE_1=/path/to/cfb-wainwright-fire-archive`. The outputs are written to this repository. To refit to the Satellogic outline:

```
cd cfb-wainwright-fire-archive-2/scripts/rubicon_ridge/model_v5
python3 v5_fit.py --scar ../../../data/rubicon-ridge/satellogic_2026-10-20/outline_for_model_v5.geojson --scar-mode union --tag satellogic
```

To rebuild model v5.1 (aircraft timeline; also reads part 1's adsb.lol trace files):

```
cd cfb-wainwright-fire-archive-2/scripts/rubicon_ridge/model_v5_1
python3 v51_air.py && python3 v51_figures.py --gif
```

Scripts that pull from services needing a key read it from an environment variable or a local file outside the repository. Never commit keys.

## Times

Alberta was on Mountain Daylight Time (UTC−6) on Oct. 17, 2026. Satellite data is in UTC. 3:29 p.m. MDT = 21:29 UTC. Alberta ACIS hourly data is in local standard time (MDT = LST + 1 h).

## Licences

| Material | Licence |
|---|---|
| Code (`scripts/`) | [MIT](LICENSE) |
| Prime Rogue's own data, analysis and write-ups, including the measurements derived from commercial imagery | [CC BY-NC 4.0](LICENSE-DATA). Credit "Prime Rogue Inc." Commercial use (including commercial media) needs written permission. |
| NASA FIRMS data | NASA open data. Credit: "NASA FIRMS (LANCE)". |
| NOAA GOES and RAP data | Public domain (U.S. government). Credit: "NOAA". |
| Copernicus Sentinel-3 data and derivatives | Free, with credit: "Contains modified Copernicus Sentinel data [2026]". |
| Environment Canada data (stations, HRDPS, radar) | Open Government Licence – Canada. Credit: "Environment and Climate Change Canada". |
| Alberta ACIS and air-quality data | Open Government Licence – Alberta. Credit: "Alberta Agriculture and Irrigation (ACIS)" / "Alberta Environment and Protected Areas". |
| PurpleAir sensor data | PurpleAir's terms of service. Credit: "PurpleAir". |
| ADS-B aircraft positions (adsb.lol open archive), including the tables derived from them in `model_v5_1_rubicon_cloud/` | [ODbL 1.0](https://opendatacommons.org/licenses/odbl/1-0/). Credit: "© adsb.lol contributors". |
| OpenStreetMap features used for reference | ODbL. Credit: "© OpenStreetMap contributors". |

**Commercial imagery is never stored here.** Prime Rogue licensed the Satellogic image of Oct. 20, 2026 through SkyFi, and the licence forbids redistribution. Only Prime Rogue's own measurements derived from it appear here, as vectors and tables with no image data. The same applies to any Planet, Vantor/Maxar or Airbus imagery.

## Contributing

Please **do not** open issues or pull requests containing documents, tips or anything that identifies you — GitHub history is permanent and public. Send material through Signal (**Rubicon49.50**) or Proton email instead; see the archive page on primerogueinc.com. Contributors are credited or protected per their preference under the Canadian Association of Journalists' ethics guidelines.

Corrections to code or methods are welcome as issues. Corrections are made in new commits, so the history shows what changed and when.
