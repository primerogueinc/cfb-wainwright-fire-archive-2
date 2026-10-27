# ADS-B positions near CFB Wainwright, Oct. 17-18, 2026 (supplied file)

`wainwright_fire_night_low_near_base_all_positions.csv` is kept here unaltered. Prime Rogue supplied it on Oct. 27, 2026. `SHA256SUMS` holds its checksum:

```
2fb3f95bfc246a45a67d3d2fae56924c4445fe14690e6ec3e47a7e3685c42602  wainwright_fire_night_low_near_base_all_positions.csv
```

## What it is

- **Size:** 1,671 aircraft positions covering 11 aircraft of 5 types (BE20, B350, C510, EC45, BK17).
- **Time span:** 15:31:28 MDT Oct. 17 to 03:57:58 MDT Oct. 18.
- **Columns:** `utc`, `mdt` (UTC−6), `hex` (ICAO 24-bit address), `registration`, `type`, `description`, `military`, `lat`, `lon`, `alt_ft` (barometric), `on_ground`, `km_from_base`.

## Provenance (checked by `scripts/rubicon_ridge/model_v5_1/v51_air.py`; results in `model_v5_1_rubicon_cloud/provenance_check.json`)

- **Source.** Every row matches a position in the [adsb.lol](https://adsb.lol) open archive (`adsblol/globe_history_2026`, readsb "prod" releases for Oct. 17 and 18). Part 1 already holds the original trace files in [`data/adsb/oct17/raw_traces/`](https://github.com/primerogueinc/cfb-wainwright-fire-archive/tree/main/data/adsb/oct17/raw_traces).
  - Every row matches an archive position with the same ICAO address to within one second.
  - Positions agree to within 0.0006°.
  - Barometric altitudes are identical in 99.9% of rows.
  - The file's times are the archive's fractional-second times rounded to whole seconds.
  - **This file adds no positions that part 1 lacked.** It is a selection from them.
- **"Base".** `km_from_base` is the distance from a round-number point at **52.80 N, 110.90 W**. A fit recovers that point, and the largest residual is 0.05 km. The point lies 4.0 km north of the v5 fire centre and 2.8 km west of Wainwright Airport (CYWV). It is not the CFB Wainwright airfield or weather station (CZWN).
- **Selection.** The file is exactly this filter applied to the archive:
  - positions within 40 km of that point;
  - numeric barometric altitude of 12,000 ft or less, which excludes positions flagged "ground";
  - a time window.
  - All 1,666 qualifying archive positions in that window are in the file. In the other five rows, the time rounded to the same second as another row, so they match an archive position that another row also matches.
  - The start cut lies between 11:20:36 MDT (the previous qualifying position) and 15:31:28 MDT, consistent with a cut at about 15:30 MDT. The end cut lies between 03:57:58 and 11:49:36 MDT on Oct. 18.
  - The 4,250 ft minimum altitude is not part of the filter. It is simply the lowest altitude the volunteer receivers heard near the base, about 2,000 ft above the ground. The "low" in the file name refers to the 12,000 ft ceiling.
- **Columns `military` and `on_ground`.** Both are `False` in every row. Military aircraft rarely broadcast ADS-B, so their absence proves nothing.

**Attribution.** The data was supplied by Prime Rogue and matches the adsb.lol open archive row for row. Data © adsb.lol contributors, [Open Database License (ODbL) 1.0](https://opendatacommons.org/licenses/odbl/1-0/).

Civil registrations and aircraft types come from public ADS-B broadcasts. Nothing in this archive says who was aboard any aircraft.
