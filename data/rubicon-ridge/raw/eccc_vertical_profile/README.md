# ECCC vertical profiles (forecast and observed), Oct. 16–18, 2026

Retrieved unaltered from Environment and Climate Change Canada's MSC Datamart on Oct. 29, 2026, about 07:45 UTC, before the Datamart's retention window drops them. Checksums are in `SHA256SUMS`.

## `forecast_RDPS_CZWN/`: forecast profiles for the CFB Wainwright airfield

ECCC's Regional Deterministic Prediction System (RDPS) forecast profiles ("ProgTephi") for station CZWN, "Wainwright_Airfield", 52.833 N, 111.100 W. Each file is one model run, with forecasts every 3 hours out to 48 hours. Variables: TT (temperature), ES (dew-point depression), HR (relative humidity), GZ (height), WD (wind direction), UV (wind speed) and WW (vertical motion), by pressure level.

| File | Model run (UTC) | Run time (MDT) | Source |
|---|---|---|---|
| `20261016_ProgTephi_00_CZWN.csv` | 00Z Oct. 16 | 6 p.m. Oct. 15 | https://dd.weather.gc.ca/20261016/WXO-DD/vertical_profile/forecast/csv/ProgTephi_00_CZWN.csv |
| `20261016_ProgTephi_12_CZWN.csv` | 12Z Oct. 16 | 6 a.m. Oct. 16 | https://dd.weather.gc.ca/20261016/WXO-DD/vertical_profile/forecast/csv/ProgTephi_12_CZWN.csv |
| `20261017_ProgTephi_00_CZWN.csv` | 00Z Oct. 17 | 6 p.m. Oct. 16 | https://dd.weather.gc.ca/20261017/WXO-DD/vertical_profile/forecast/csv/ProgTephi_00_CZWN.csv |
| `20261017_ProgTephi_12_CZWN.csv` | 12Z Oct. 17 | 6 a.m. Oct. 17 | https://dd.weather.gc.ca/20261017/WXO-DD/vertical_profile/forecast/csv/ProgTephi_12_CZWN.csv |

Run time is the model's start time. ECCC publishes the output a few hours later.

The files don't state a unit for wind speed (UV). Prime Rogue reads it as knots, which is consistent with the observed soundings; confirmation from ECCC is pending.

## `observed_soundings/`: weather-balloon soundings

ECCC observed upper-air soundings ("ObsTephi") at 00Z and 12Z (6 p.m. and 6 a.m. MDT), Oct. 16–18. Stations:

- CWSE: Stony Plain, near Edmonton (the closest to Wainwright);
- CYQD: The Pas, Man.;
- CWVK: Vernon, B.C.

Source: `https://dd.weather.gc.ca/YYYYMMDD/WXO-DD/vertical_profile/observation/csv/ObsTephi_HH_STATION.csv`

Data: Environment and Climate Change Canada, Open Government Licence – Canada.
