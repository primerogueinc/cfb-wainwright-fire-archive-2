# Polar-orbiter passes (VIIRS / MODIS), CFB Wainwright fire, Oct. 17, 2026, 13:00-18:00 MDT

Source for pass geometry: SGP4 propagation of current Space-Track TLEs (via SatNOGS DB, fetched Oct. 27,
checksums verified, 5-15 days back-propagated) for Suomi-NPP, NOAA-20, NOAA-21, Aqua and Terra.
Cloud context comes from the GOES-18 Full Disk point series already in `raw/goes/`.
Script: `scripts/rubicon_ridge/polar_visible.py` (`passes`, then `extract` and `render`).

## Bottom line

**No polar-orbiter imagery could be examined, so this folder makes no plume claim, for or against.**
The VIIRS buckets (`noaa-nesdis-snpp-pds`, `noaa-nesdis-n20-pds`, `noaa-nesdis-n21-pds`) and NASA
GIBS/Worldview, CMR and LAADS were all refused by this session's network policy (HTTP 403). No
figures were rendered. What the pass geometry does establish:

- **No polar pass imaged the fire during its intense phase (about 15:30-15:50 MDT).** The only pass near
  that window, NOAA-21 at 15:35 MDT, misses: the fire lies about 25 km beyond the west edge of the swath
  (scan angle 56.3° vs the 56.06° VIIRS limit). Even if the swath reached it, the smoke was blowing SE, away
  from the swath edge. Confirm from the granule's GEO file when access exists.
- **After the front, the only pass is Aqua MODIS at 16:55 MDT** (36° scan). GOES-18 shows cold cloud over
  the fire at that time (10.3 µm 260 K) and no fire heat after 15:50, so a plume is unlikely to be visible.
  The sun is low (12°).
- **The best lead is Aqua MODIS at 15:18 MDT**, about 10 minutes before the front. It is the only in-swath
  pass with GOES-18 indicating a mostly clear view at the fire (10.3 µm 285 K, FDC mask 100, "processed, no fire", meaning it was not cloud-flagged).
  GOES-18 first detected the fire at 15:10 (76 MW). The view is near the swath edge (53° scan), so 250 m
  pixels are about 1.0 × 0.5 km and 500 m true-colour pixels about 2 × 1 km. Any plume then would be the
  pre-front one, carried roughly ENE on the WSW wind (RAP 258° at 15:00), not the post-front SE plume.
  Not checked: whether Aqua MODIS was still producing data on Oct. 17, 2026. FIRMS returned no MODIS
  detections that day, which could be cloud, a small fire, an edge-of-swath miss or no data.
- **NOAA-20 at 14:50 MDT is the second lead.** VIIRS detected the fire (8 I-band pixels, FRP 6-12 MW),
  but GOES-18 shows partly cloudy or thin high cloud at the fire (10.3 µm 257 K, rising from about 230 K).
  This was a small, pre-front fire. Any plume would be faint and pointing roughly ENE.
- **The earlier passes are clouded out.** NOAA-20 at 13:10, NOAA-21 at 13:54 (near nadir) and S-NPP at 14:31
  all coincide with cold cloud at the fire on GOES-18 (10.3 µm 225-233 K). This is consistent with no
  FIRMS detections from NOAA-21 or S-NPP that day, even though Landsat 8 saw the fire at 12:17.
- **Sentinel-3 OLCI/SLSTR and Landsat have no pass in the window.** Their daytime passes were 11:30
  (S3B), 12:09 (S3A) and 12:17 MDT (Landsat 8); see `raw/sentinel3_frp/README.md`.

## Pass list

All passes over the fire whose closest approach falls between 12:00 and 18:40 MDT. "Swath" means the fire
is inside the instrument's scan limit (VIIRS ±56.06°, MODIS ±55°). Pixel sizes are for the finest band:
VIIRS I-band 375 m, MODIS 250 m. They are given as along-scan × along-track and are approximate. The
GOES-18 10.3 µm values are for the fire pixel, from the scan nearest the pass. Values of about 225-240 K
mean cold, high cloud; values near 285 K mean clear ground before the front.

| MDT | Satellite | Scan / view zenith | Swath | Pixel (km) | Sun elev. | G18 10.3 µm at fire | Status |
|---|---|---|---|---|---|---|---|
| 12:14 | NOAA-21 | 57.1° / 71.8° | no | – | 27° | 234 K | before window |
| 12:51 | S-NPP | 53.2° / 64.8° | yes | 0.59 × 0.72 | 28° | 229 K | before window, cloud |
| 13:10 | NOAA-20 | 45.4° / 53.7° | yes | 0.34 × 0.57 | 28° | 229 K | cloud |
| 13:31 | Terra | 63.9° / 85.5° | no | – | 28° | 233 K | out of swath |
| 13:42 | Aqua | 63.8° / 84.6° | no | – | 28° | 229 K | out of swath |
| 13:54 | NOAA-21 | 9.3° / 10.6° | yes | 0.40 × 0.38 | 27° | 225 K | cloud; no FIRMS detection |
| 14:31 | S-NPP | 17.1° / 19.4° | yes | 0.43 × 0.39 | 26° | 233 K | cloud; no FIRMS detection |
| **14:50** | **NOAA-20** | 38.5° / 44.8° | yes | 0.49 × 0.50 | 24° | 257 K | **lead**: fire detected (FIRMS 14:49), thin/partial cloud |
| **15:18** | **Aqua (MODIS)** | 53.0° / 62.4° | yes (edge) | 1.01 × 0.47 | 22° | 285 K | **best lead**: mostly clear at fire, pre-front |
| 15:35 | NOAA-21 | 56.3° / 70.2° | no (≈25 km outside) | – | 21° | 278 K | misses intense phase |
| 16:12 | S-NPP | 59.6° / 77.1° | no | – | 17° | 268 K | out of swath |
| 16:32 | NOAA-20 | 61.3° / 82.9° | no | – | 15° | 260 K | out of swath |
| 16:55 | Aqua (MODIS) | 35.9° / 40.6° | yes | 0.42 × 0.32 | 12° | 260 K | post-front cloud |
| 18:15 | Terra | 64.1° / 87.1° | no | – | 1° | 278 K | out of swath, dusk |

Full details, including UTC times, satellite azimuth, sub-satellite point, ground distance from the track
and TLE age, are in `passes.csv`.

## Limits

- **Times.** Closest-approach times are good to about ±1 min, given the TLE age. VIIRS and MODIS image the
  fire's scan line within a few seconds of closest approach. Check: the NOAA-20 pass is predicted at 20:50:32
  UTC, and FIRMS lists 20:49 (its granule time), with a reported pixel size of 0.48 × 0.48 km against the
  predicted 0.49 × 0.50.
- **Swath edge.** The NOAA-21 15:35 miss (56.3° vs 56.06°, about 25 km outside) is larger than the expected
  cross-track error, but it is a computed result, not one confirmed against a granule.
- **Cloud.** These values come from GOES-18 at a 65° view angle. High cloud seen over the fire on GOES-18 is
  physically about 20 km SW of the fire (parallax at about 10 km height). However, the burn-box median and the
  background 25 km WSW show the same values, so the cloud signal is regional. Thin cirrus at about 230 K can
  still let a polar sensor see some ground. The "cloud" rows are therefore likely, not certain.
- **No imagery examined.** Every statement above about plume visibility is an expectation, not an
  observation.

## To finish (once the hosts are reachable)

1. Allow `noaa-nesdis-n20-pds.s3.amazonaws.com` (and the snpp/n21 buckets). Then run
   `python3 -I scripts/rubicon_ridge/polar_visible.py extract <this dir> NOAA-20 S-NPP NOAA-21`, then `render`.
   The code is written but has not been run against real granules. It reads only a ~220 km block of each
   granule over HTTPS range reads and writes small `crops_*.npz` files. It stores no raw granules.
2. For Aqua MODIS 15:18 and 16:55, the quickest route is Worldview: layer
   `MODIS_Aqua_CorrectedReflectance_TrueColor`, date 2026-10-17. LAADS MYD02QKM/MYD02HKM + MYD03 would be
   needed for a figure in this style. This script does not cover MODIS.

Credit: NOAA/NASA VIIRS via NOAA Open Data on AWS (when rendered); pass geometry from Space-Track TLEs.
