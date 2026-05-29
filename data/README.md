# `data/`

All data for the project. Top-level layout:

| Subdirectory | What it holds |
| --- | --- |
| [`raw_data/`](raw_data) | Acquisition time series, one CSV per campaign day, plus the extraction windows used by the pipeline. |
| [`derived_data/`](derived_data) | Products computed from the raw data: per-isenthalp $p$-$T$ pairs, Joule-Thomson coefficients, and capillary-sizing design calculations. |
| [`cernox_calibration_data/`](cernox_calibration_data) | Lake Shore calibration files for the Cernox temperature sensors. |
| [`literature/`](literature) | Historical Joule-Thomson measurements digitised for the comparison figure. |
| [`capillary_design/`](capillary_design) | Pressure-drop design calculation for the abandoned 10 m capillary geometry. |
| [`metadata/`](metadata) | Descriptive metadata: column dictionary, file manifest, sensor provenance, campaign log, DAQ stack. |
| `all_measurements.csv` | All raw campaign CSVs concatenated into a single file for convenience. |
