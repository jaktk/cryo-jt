# `p_T_pairs`

This directory contains the steady-state pressure-temperature(-composition)
points extracted from the raw time series in [`../../raw_data`](../../raw_data)
by [`get_pT_pairs.py`](../../../src/data_analysis/get_pT_pairs.py), one CSV
per isenthalp.

Each row is the mean over one steady-state window of
[`extraction_windows.json`](../../raw_data/extraction_windows.json) (5 samples
at 1 s). For every raw column `X` the files also hold its standard deviation
`X_STD` and the expanded uncertainty of the mean `X_EXP_UNC` = 2 `X_STD`/√N
(k = 2, about 95 % level of confidence).

File names follow `<fluid>_<T_in>K_<p_in>MPa.csv`, where `T_in` and `p_in` are
the mean inlet temperature and pressure of the isenthalp, rounded to K and
MPa. The isenthalps marked `"excluded"` in `extraction_windows.json` are not
used in the published analysis and are not extracted.

[`index.json`](index.json) is the manifest read by
`calculate_jt_coefficients.py` (file name, fluid, mean inlet temperature). It
is not written by `get_pT_pairs.py`, so it must be updated by hand when an
isenthalp is added or removed.
