# `data/`

All data for the project. Top-level layout:

| Subdirectory | What it holds |
| --- | --- |
| [`raw_data/`](raw_data) | Acquisition time series, one CSV per campaign day, plus the extraction windows used by the pipeline. |
| [`derived_data/`](derived_data) | Products computed from the raw data: per-isenthalp $p$-$T$ pairs, Joule-Thomson coefficients with uncertainties, correction estimates, generated LaTeX tables, and capillary mass-flow calculations. |
| [`cernox_calibration_data/`](cernox_calibration_data) | Lake Shore calibration files for the Cernox temperature sensors. |
| [`literature/`](literature) | Historical Joule-Thomson measurements digitised for the comparison figure. |
| [`capillary_design/`](capillary_design) | Pressure-drop design calculation for the abandoned 10 m capillary geometry. |
| [`metadata/`](metadata) | Descriptive metadata: column dictionary, file manifest, sensor provenance, campaign log, DAQ stack. |
| `all_measurements.csv` | One raw sample per steady-state measurement point (126 points over 27 isenthalps, numbered in column `no`), compiled as an overview. It is not a concatenation of the raw time series and is not read by the pipeline; its isenthalp numbering and time stamps do not all match `raw_data/extraction_windows.json`. |
