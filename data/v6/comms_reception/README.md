# v6 working material: communications, reception and weather alerts (Oct. 29, 2026)

Working material for the CFB Wainwright fire (Oct. 17, 2026). It is **unpublished v6 analysis**: models and summaries, not findings of fact. Analysis, maps and text: Prime Rogue Inc., CC BY-NC 4.0.

| Folder | What's in it |
|---|---|
| `reports/` | PDFs. The full technical radio and phone reach report (`Prime_Rogue_Radio_and_phone_reach_CFB_Wainwright.pdf`), the phone-reception model, the WiGLE review, the licensed-radio comparison across six Army training areas, the Oct. 16–17 weather-alerts brief, and the knots vs km/h brief |
| `ised_licences/` | Civilian licensed stations near the burn (ISED Spectrum Management System extracts, Oct. 5, 2026); comparison across training areas; public map |
| `wigle/` | Summary and map of WiGLE cell and Wi-Fi sightings. Raw WiGLE records are not republished here, and no network names are recorded |
| `cell_model/` | Phone-reception model: tower list, scripts, results (JSON), figures |
| `radio_model/` | Radio reach model (handheld, vehicle, base station, Wainwright fire VHF): scripts, sensitivity runs, weather and active-fire scenarios, figures |
| `weather_alerts/` | Every ECCC and aviation weather message for Alberta, Oct. 16–17 (226 events), parsed CAP alert flags, scripts and charts |

## Key assumptions (see the technical report)
- **Radios are generic, not DND's equipment.** Low VHF 50 MHz handheld 5 W / vehicle 50 W, receiver −113 dBm, 8 dB location variability. The garrison base-station site (52.825 N, 110.885 W, 20 m mast) is assumed.
- **Terrain:** Copernicus GLO-30 surface model, downloaded by the scripts. **Refraction:** standard k = 4/3, and k = 1.2 from the ECCC Stony Plain sounding at 6 p.m. Oct. 17. **Active-fire case:** an assumed +6 dB loss inside the burn.
- **Pairs-range location** (52.7695 N, 110.9005 W) is our best estimate, ±300 m.

## Deliberately left out
- Frequencies and exact positions of DND's licensed microwave links (public in ISED data, but not needed here).
- Raw WiGLE records (WiGLE terms) and any Wi-Fi network names.

## Reproducing
The scripts were run from a working directory with sibling folders (`../v6/s2_lens/search.json` is a Sentinel-2 Earth Search query result; `../radio/` holds the ISED station table). Adjust the paths at the top of each script. Python with numpy, scipy, rasterio, geopandas, pyproj and matplotlib.

## Sources
ISED Spectrum Management System extracts (Open Government Licence – Canada); Copernicus DEM GLO-30 and Sentinel-2 (contains modified Copernicus data, 2026); ECCC MSC Datamart bulletins, CAP alerts and soundings (OGL – Canada); NAV CANADA/ECCC aviation bulletins via the Datamart; WiGLE.net (summary only).
