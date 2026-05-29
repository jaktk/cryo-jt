# `data/metadata/`

Descriptive metadata for the dataset: what each file is, what every column
means, which sensors were used, when the data were acquired, and on what
hardware/software. The files are static descriptions — they document the
dataset but are not consumed by the analysis pipeline.

Operational inputs to the pipeline (extraction windows and the index of
derived isenthalp files) live next to the data they act on:

- `data/raw_data/extraction_windows.json` — time stamps used by
  [`get_pT_pairs.py`](../../src/data_analysis/get_pT_pairs.py) to cut each
  steady-state measurement window out of the raw time series.
- `data/derived_data/p_T_pairs/index.json` — manifest of the per-isenthalp
  CSVs, used by
  [`calculate_jt_coefficients.py`](../../src/data_analysis/calculate_jt_coefficients.py)
  to look up the fluid for each file.

## Contents

| File | What it is |
| --- | --- |
| [`README.md`](README.md) | This file. |
| [`codebook.csv`](codebook.csv) | Variable dictionary: per-column type, allowed values, unit, and description for every CSV produced or consumed by the pipeline. |
| [`file_list.csv`](file_list.csv) | Manifest of the data files in the repository, with a one-line description and license for each. |
| [`sensors.json`](sensors.json) | P&ID tag ↔ sensor model ↔ serial number ↔ calibration provenance for the principal measurement chain (TT101, TT102, PT101, PT102, GA016). |
| [`campaigns.csv`](campaigns.csv) | Measurement campaign log: date, fluid or mixture, mean composition or purity, raw-data filename. |
| [`daq.json`](daq.json) | Data-acquisition stack: SCADA software, temperature acquisition module, transport. |

## How to use it

- Looking at an unfamiliar column in any CSV? Find it in `codebook.csv`.
- Want to know which sensor produced `TT101/K` and how its calibration
  uncertainty was derived? See `sensors.json`.
- Want to know which raw file corresponds to a given campaign date or
  fluid? See `campaigns.csv`.
- Reproducing the acquisition setup? See `daq.json`.
