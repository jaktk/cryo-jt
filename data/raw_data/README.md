# `data/raw_data/`

Acquisition time series, sampled at 1 s. One CSV per measurement campaign
day, named `YYYY-MM-DD.csv`. Column definitions are in
[`../metadata/codebook.csv`](../metadata/codebook.csv); the fluid or
mixture measured on each day is listed in
[`../metadata/campaigns.csv`](../metadata/campaigns.csv).

`extraction_windows.json` lists the steady-state start/end time stamps used by
[`get_pT_pairs.py`](../../src/data_analysis/get_pT_pairs.py) to slice each
isenthalp's measurement points out of the raw time series.
