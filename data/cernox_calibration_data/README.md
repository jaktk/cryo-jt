# `data/cernox_calibration_data/`

Lake Shore Cryotronics calibration files for the Cernox resistive
temperature sensors used as TT101 and TT102. One subdirectory per sensor
serial number.

| Serial | Used as | Lake Shore model | Calibration range |
| --- | --- | --- | --- |
| [`X93303`](X93303) | TT102 (outlet, in fluid contact) | CX-1050-SD-HT-1.4L | 1.4 - 325 K |
| [`X115143`](X115143) | TT101 (inlet, on the T-piece) | CX-1080-CU-HT-20L | 20 - 325 K |

The TT101 / TT102 ↔ serial mapping is also recorded in
[`../metadata/sensors.json`](../metadata/sensors.json).

Per-sensor file types (one of each in every subdirectory):

| Extension | Content |
| --- | --- |
| `*.cof` | Chebyshev polynomial coefficients $T(R)$ on the calibrated range |
| `*.dat` | Calibration point list (resistance and temperature pairs) |
| `*.tbl` | Calibration look-up table |
| `*.curve` | Lake Shore-format curve file for direct upload to a 218 / 224 / 240 temperature monitor |
| `*.234`, `*.330`, `*.340`, `*.34A`, `*.91C` | Lake Shore-format curve files for the corresponding instrument models |
| `*.pdf` | Scanned calibration certificate |
