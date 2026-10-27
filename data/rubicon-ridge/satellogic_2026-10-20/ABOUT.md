# Satellogic 0.7 m image, Oct. 20, 2026: Oct. 17 burn, CFB Wainwright

Derived vectors and tables only. The image is licensed (Satellogic via SkyFi) and is **not** in this repository: no rasters, crops or chips of it are stored here. Credit any figure made from it as "Imagery: Satellogic via SkyFi, Oct. 20, 2026".

## Image

- Satellogic NewSat-49, 2026-10-20 18:56:43 UTC (12:56 MDT), three days after the fire.
- 0.7 m GSD. Four bands (red, green, blue, NIR) as uint16 TOA reflectance x 1e4. EPSG:32612, 8928 x 8928 px.
- Off-nadir angle 24.2°, satellite azimuth 36°, sun elevation 27°, sun azimuth 179°, so shadows fall to the north.
- Footprint: about 503,580–509,830 m E and 5,843,360–5,849,610 m N (`image_footprint.geojson`). Range 11 is outside it, about 1 km west of the image edge. The south-east tip of the S2/Landsat outline is also outside it: 24 ha of that outline's 1,079 ha.
- The image has no SWIR band, so no NBR/dNBR is possible. All burn mapping uses visible and NIR only.

## 1. Georeference

The image was compared with Sentinel-2 L2A from Oct. 19, 2026, using bands B02, B03, B04 and B08. Each band was reprojected to a 10 m grid and phase-correlated over 24 tiles of 1.6 km.

- Result: features in the Satellogic image lie **1.7 m west and 0.3 m north** of the same features in Sentinel-2. Tile MAD is 0.25–0.5 m. All four bands agree.
- Sentinel-2 from Oct. 4 agrees (west 0.5–2 m). So does Landsat from Oct. 24 (west 3.8 m, south 1.5 m, at 30 m, MAD 1.5 m).
- OSM road centrelines give 2.0 m west. The OSM north component is poorly constrained because most of the roads run north-south.
- **Correction applied:** every vector in this folder is shifted +1.7 m E and −0.3 m N into the Sentinel-2 frame, which is the frame of `scar_patches_revised` and the model.
- The offset is smaller than Sentinel-2's own absolute geolocation uncertainty, so absolute positions carry roughly ±3–5 m. Details are in `georeference_offset.json`.

## 2. Burn outline

### Method

1. **Training labels.** Burned pixels were taken from at least 60 m inside the S2/Landsat "main Oct 17 run". Unburned pixels were taken from at least 150 m outside every scar patch and outside the earlier 2026 burns (Sept. 18 PDA 2, Sept. 20, Oct. 14).
2. **Classifier.** A gradient-boosting classifier was run at full resolution. Its inputs were NIR, red, visible mean, NDVI and (R−B)/I, plus 7.7 m and 21.7 m local means and 7.7 m texture.
   - Holdout recall when training on one half of the image (east or west) and testing on the other: burned 0.96, unburned 0.90–0.97.
   - The texture and context features were added because an earlier quadratic classifier marked leafless aspen bluffs and shaded dune slopes as burned.
3. **Cleaning.** Pixels with P(burned) ≥ 0.5 were kept. Then:
   - a 3.5 m majority filter was applied;
   - earlier burns (buffered 30 m) and OSM water were removed;
   - only components connected to the main-run core were kept, with gaps up to about 12 m bridged;
   - burned specks under 25 m² were dropped;
   - unburned holes of 0.1 ha or more were kept as islands.
4. **Excluded burns.** The Oct. 14, Sept. 20 and Sept. 18 burns look spectrally the same as Oct. 17 ground three to six days later. They are excluded by their pre-Oct. 17 Sentinel-2 outlines, not by spectrum. The April 20, 2026 planned burn has six months of regrowth and separates clearly by NIR (NDVI 0.235 vs 0.13 on fresh burn). Where the Oct. 17 fire re-burned it, the ground is classed as burned.

### Results

`burn_outline_oct17_satellogic.geojson`:

| | Area |
|---|---|
| Net burned (excluding 52 unburned or water islands ≥ 0.1 ha) | **1,049.8 ha** |
| Envelope including islands | 1,085.2 ha |
| Islands | 35.4 ha |
| Sensitivity, P ≥ 0.35 | 1,063 ha net |
| Sensitivity, P ≥ 0.65 | 1,037 ha net |

The sensitivity outlines are in `burn_outline_sensitivity.geojson`.

### Comparison with `data/infrastructure/scar_patches_revised.geojson` (main Oct 17 run)

The S2/Landsat outline has no islands, so the comparison uses the envelope. Inside the image footprint:

- S2/Landsat outline: 1,055.1 ha. Satellogic envelope: 1,085.2 ha, which is 30.1 ha (2.9%) larger.
- Burned in Satellogic only: 42.2 ha. Inside the S2/Landsat outline only: 12.1 ha. IoU 0.95.
- Polygons ≥ 0.25 ha are in `outline_change_vs_scar_patches_revised.geojson`.

Edge displacement was sampled every 20 m along the old edge (`edge_displacement_vs_scar_patches_revised.csv`; + means the new edge lies outside the old one):

- Median +3.7 m; median absolute 9.5 m.
- 52% of points within 10 m, 77% within 20 m.
- 90th percentile of absolute displacement: 50 m.

Where the edges moved:

- **West and south-west lake shores:** moved out by a median of 78 m in the SW sector, up to about 270 m. The largest single gain is 22.3 ha, burned grass between the two western saline lakes (centred about 505,570 E / 5,845,370 N) that the 10–30 m outline cut across.
- **North tip and north:** a few metres to tens of metres. The new outline follows the dune-slope edge more tightly and adds about 1.7 ha between the tracks just south of the wet depression at the north tip.
- **East and north-east:** median within ±1 m, up to −42 m where the old outline took in unburned aspen pockets.
- **South-east (east of the boundary road):** about 1.8 ha inside the old outline near the image edge is mapped unburned, and aspen stands of up to 2.5 ha appear unburned or crown-unburned and are mapped as islands. About 1.9 ha of burned ground is added nearby.
- **The former "inside" of the old outline:** the 52 islands are mostly wet or dry ponds (the largest is a 9.2 ha pond at about 504,530 E / 5,846,390 N), unburned aspen clumps and grass pockets.

`outline_for_model_v5.geojson` is the same net outline (EPSG:4326, S2 frame) for the v5 refit. Note that it stops at the image edge in the south-east.

## 3. Edges

See `edge_findings.csv`, `edge_runs_along_features.geojson` (outline runs within 40 m of a road, track or the OSM PDA 2 boundary) and figure F3.

Measurement method:

- Run lengths are the length of outline lying within 40 m of a feature.
- Offsets come from median image cross-profiles across the feature. Positive offsets point toward the burn.
- A "burn-out" was searched for as a narrow band, within 15 m of the edge on the burned side, that is at least 8% darker in NIR than the burned interior.

Findings:

- **E1 – north of the road head.** About 740 m of outline runs parallel to the N-S road. The burn stops 27–29 m west of the road centreline, at the west margin of an unburned 20–25 m corridor (NDVI about 0.24, like unburned grass). The OSM-mapped PDA 2 boundary lies a median of about 7 m from the outline, with about 1.1 km within 40 m.
- **E1b – north tip.** The burn stops at the N-S track. A very dark, wet depression immediately west of the track burned to black ash (NDVI about 0.08, NIR well below the burn median).
- **E2 – Range 2 access road.** About 470 m of outline is within about 9 m of the road centreline, on the road's north side. South of the road an unburned grass lens about 400 m long separates it from the main burn.
- **E3 – road head southward.** The burn edge is on the west shoulder of the N-S road for about 90 m south of the junction. From about 90 m south of the junction, burned ground lies on both sides of the road, so the fire crossed it.
- **E4 – east edge.** About 1.4 km long, nearly straight (bearing 178°). It sits at the west margin of the base-boundary road allowance, with about 1.47 km within 40 m of the OSM PDA 2 boundary. The travelled track is 35–75 m further east. Nothing burned east of it north of the farm track (about 52.752 N).
- **E5 – south of a farmyard.** Over about 0.7–0.9 km, the burn stops at an E-W farm track. South of about 52.752 N, burned ground lies on both sides of the boundary road.
- **E6 – south-east gravel road dog-leg.**
  - SE arm: the burn stops at the road's north-east side.
  - NE arm: straight for about 320 m (RMS 2 m). The 20–30 m road strip is unburned but burned ground lies on both sides.
- **E7 – north-west flank and E8 – north-east flank.** Long flanks (3.6 km and 3.0 km of outline) with no track. Neither is straight in an engineered sense (RMS 77 m and 42 m from a line).
- **Burn-outs.** No continuous straight black line along any track was found. One short NIR dip (about 3.5 m) along the E6 SE arm, and one on an OSM track at the south shore, are inconclusive.

## 4. Points

`points_range2_range15_roadhead.csv` uses the official grid refs as given. "Outline" metrics use the net outline polygon; "pixel" metrics use the 0.7 m burned-ground mask. Both agree to within 1 point.

| Point | At point | % burned within 50 m | % burned within 150 m | Distance to edge | Nearest burn |
|---|---|---|---|---|---|
| Range 2 | unburned | 20 | 31 | 16 m | due north, across the access road |
| Range 15 | unburned | 27 | 44 | 12 m | WSW, across the N-S road |
| Road head | unburned | 35 | 43 | 6 m | SW, the west shoulder |
| Range 11 | outside image | – | – | – | – |

Ground descriptions:

- **Range 2** is on the south shoulder of its access road, about 80–90 m west of the N-S road junction, in open grass beside a pale graded or gravel area.
- **Range 15 and the road head** are on open grass/gravel at the junction of the N-S road and the E-W track, east side.

## 5. Surface classes (qualitative)

`surface_classes_qualitative.geojson` covers the burned ground on a 2.8 m grid, with a 14 m majority filter and a 0.1 ha minimum:

| Class | Area | Share |
|---|---|---|
| Black ash | 200 ha | 20% |
| Intermediate dark-grey ash | 558 ha | 55% |
| Grey/light ash | 213 ha | 21% |
| Standing vegetation (scorched or crown-unburned) | 46 ha | 4% |
| Bright track/road/bare surface | 1 ha | <1% |

Caveats:

- The rules are fixed TOA thresholds; see the script note below and `summary_stats.json`.
- These are appearance classes, not calibrated burn severity.
- True white ash is rare or unresolvable: bright, low-NDVI pixels make up under 0.2% of the burn. "Grey/light ash" means lighter-toned, low-vegetation burned ground and may include exposed sand.
- With the sun at 27°, north-facing dune slopes are darker and bias toward "black ash".
- Leafless aspen over burned ground cannot be told apart from scorched standing stems.

## 6. Objects

`objects_near_ranges_and_tracks.csv` and `.geojson` cover the AOI: within 300 m of Range 2, Range 15 and the road head, and within 25 m of roads and tracks in or near the burn.

- **Method.** Candidates came from brightness and colour-fringe anomalies against a 17.5 m background, followed by full-resolution visual checks. Only objects are recorded.
- **Count.** 16 objects:
  - 2 vehicle-sized objects on the N-S road (O1, O2) whose band-to-band colour smear is consistent with movement at acquisition;
  - 2 small objects on the west shoulder at the road head (O6, O7);
  - 2 large bright objects (O4 at the gravel pad at the west end of the Range 2 access road; O8 about 85 m SE of Range 15). Each is up to 10–14 m (saturated) with a shadow on its north side: a structure, trailer or large vehicle, not resolvable;
  - 1 bright vehicle- or trailer-sized object on the N track (O15);
  - 7 small items (O3, O5, O9, O14, and the pairs O10/O11 and O12/O13);
  - 1 low-confidence possible object by the east road (O16).
- **Location.** None is inside the burn.
- **Not found.** No equipment, vehicles or structures were found on the tracks inside the burn.
- **Limits.** At 0.7 m, objects under about 2 m, dark objects in shade and objects under trees can be missed. Bright objects saturate, so their sizes are upper bounds. The three-day gap means the objects describe Oct. 20, not Oct. 17.

## Files

| File | Content |
|---|---|
| `burn_outline_oct17_satellogic.geojson` | Final net outline with islands as holes (EPSG:4326, S2 frame) |
| `outline_for_model_v5.geojson` | The same geometry, packaged for the v5 refit |
| `burn_outline_sensitivity.geojson` | Outlines at P ≥ 0.35 and P ≥ 0.65 |
| `unburned_islands.geojson` | 52 islands of 0.1 ha or more |
| `outline_change_vs_scar_patches_revised.geojson`, `edge_displacement_vs_scar_patches_revised.csv` | Comparison with the S2/Landsat outline |
| `edge_findings.csv`, `edge_runs_along_features.geojson` | Edge analysis |
| `points_range2_range15_roadhead.csv` | Point results |
| `surface_classes_qualitative.geojson` | Qualitative surface classes |
| `objects_near_ranges_and_tracks.csv`, `.geojson` | Objects |
| `georeference_offset.json`, `summary_stats.json`, `image_footprint.geojson` | Georeference, statistics and footprint |

The processing scripts (s04e classifier, s05 outline, s07–s15) were run from the analysis scratchpad against the licensed image. They are not in the repo because they read the licensed files directly. Figures for internal review (F1–F4) are likewise kept outside the repo.

## General caveats

- The image is from three days after the fire. Any change between Oct. 17 and Oct. 20 is included, such as suppression or mop-up disturbance and wind-blown ash.
- The 24° off-nadir view displaces tall objects (trees, structures) by a few metres toward the south-west. This does not matter for ground-level burn edges.
- Wet, dark ground and deep shadow can mimic black ash. The connectivity rule limits but does not eliminate this.
- OSM features (ODbL) are traced from aerial imagery of unknown accuracy, typically 1–5 m, but up to tens of metres for some tracks.
