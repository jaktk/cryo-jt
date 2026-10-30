# `data/raw_data/`

Acquisition time series, sampled at 1 s. One CSV per measurement campaign
day, named `YYYY-MM-DD.csv`. Column definitions are in
[`../metadata/codebook.csv`](../metadata/codebook.csv); the fluid or
mixture measured on each day is listed in
[`../metadata/campaigns.csv`](../metadata/campaigns.csv).

`extraction_windows.json` lists the steady-state start/end time stamps used by
[`get_pT_pairs.py`](../../src/data_analysis/get_pT_pairs.py) to slice each
isenthalp's measurement points out of the raw time series.

`2021-03-04.csv` covers 11:10:00 to 16:06:29: the published nitrogen
isenthalp at 150 K/6 MPa measured with decreasing outlet pressure, followed by
the return run with increasing pressure used for the hysteresis check
(`hysteresis_check` stamps in `extraction_windows.json`). Up to version 1.3.1
the file ended at 14:20:00, before the return run.
