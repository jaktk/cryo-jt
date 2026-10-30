"""
In-situ purity of the pure fluids from the gas analyser GA016 (BGA244HP in
purity mode), for the Supporting Information (Section S1).

For every pure-fluid campaign day in data/metadata/campaigns.csv the
`impurity` column of the raw time series is filtered:
  * null frames (impurity exactly zero: no valid reading) are dropped;
  * frozen frames are dropped: runs of FROZEN_RUN or more consecutive identical
    readings. At 1 Hz a working analyser never repeats a value for that long
    (the longest run on the days with a live signal is about 130 frames), while
    on some days the output stays fixed at one value for hours.
The day-mean purity is the mean of 1 - impurity over the remaining frames;
the per-gas value pools the frames of all days of that gas, with the expanded
uncertainty taken as twice the standard deviation of the per-frame purity.

Outputs (data/derived_data/purity):
  * purity_per_day.csv   frames recorded, null, frozen, and valid, day mean
  * purity_per_gas.csv   pooled purity, frames, days, and U (k = 2)
  * si_purity_tables.tex the two tables of SI Section S1

No REFPROP needed.
"""
import os
import numpy as np
import pandas as pd
from repo_root import REPO_ROOT

FROZEN_RUN = 300   # frames (5 min at 1 Hz)
OUT_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'purity')
GAS_NAMES = {'N2': 'nitrogen', 'Ar': 'argon', 'He': 'helium-4', 'Ne': 'neon'}


def frozen_mask(values, min_run=FROZEN_RUN):
    """ True for frames inside a run of at least min_run identical values. """
    mask = np.zeros(len(values), dtype=bool)
    start = 0
    for i in range(1, len(values) + 1):
        if i == len(values) or values[i] != values[start]:
            if i - start >= min_run:
                mask[start:i] = True
            start = i
    return mask


def main():
    campaigns = pd.read_csv(os.path.join(REPO_ROOT, 'data', 'metadata', 'campaigns.csv'))
    campaigns = campaigns[~campaigns['is_mixture']]
    os.makedirs(OUT_DIR, exist_ok=True)

    days, frames = [], []
    for _, c in campaigns.iterrows():
        raw = pd.read_csv(os.path.join(REPO_ROOT, 'data', c['raw_data_file']), usecols=['impurity'])
        imp = raw['impurity'].to_numpy(dtype=float)
        null = imp == 0
        frozen = frozen_mask(imp) & ~null
        valid = ~null & ~frozen
        purity = 1.0 - imp[valid]
        days.append({'date': c['date'], 'gas': GAS_NAMES[c['fluid']], 'frames': len(imp),
                     'null': int(null.sum()), 'frozen': int(frozen.sum()), 'valid': int(valid.sum()),
                     'mean_purity': purity.mean() if valid.any() else np.nan})
        frames.append(pd.DataFrame({'gas': GAS_NAMES[c['fluid']], 'purity': purity}))
    days = pd.DataFrame(days)
    frames = pd.concat(frames, ignore_index=True)

    gases = []
    for gas, g in frames.groupby('gas', sort=False):
        d = days[(days['gas'] == gas) & (days['valid'] > 0)]
        gases.append({'gas': gas, 'purity': g['purity'].mean(), 'frames': len(g),
                      'days': len(d), 'U_k2': 2 * g['purity'].std()})
    gases = pd.DataFrame(gases)

    days.to_csv(os.path.join(OUT_DIR, 'purity_per_day.csv'), index=False)
    gases.to_csv(os.path.join(OUT_DIR, 'purity_per_gas.csv'), index=False)

    def n(x):
        return f"{x:,}".replace(',', r'\,')

    tex = [r"\begin{table}[h]", r"  \centering",
           r"  \caption{\label{tab:si_purity_day}Per-campaign in-situ purity of the pure fluids measured "
           r"with the GA016 binary gas analyzer (purity mode): frames recorded, frames discarded as null "
           rf"(zero impurity) or frozen (runs of at least {FROZEN_RUN} identical readings), valid frames $N$, "
           r"and the mean of $1-\text{impurity}$ over the valid frames.}",
           r"  \begin{tabular}{llrrrrc}", r"    \toprule",
           r"    date & gas & recorded & null & frozen & $N$ & day-mean purity \\", r"    \midrule"]
    for r in days.itertuples():
        mean = "--" if np.isnan(r.mean_purity) else f"{r.mean_purity:.4f}"
        tex.append(f"    {r.date} & {r.gas} & {n(r.frames)} & {n(r.null)} & {n(r.frozen)} & "
                   f"{n(r.valid)} & {mean} \\\\")
    tex += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}", "",
            r"\begin{table}[h]", r"  \centering",
            r"  \caption{\label{tab:si_purity_gas}Pooled in-situ mole-fraction purity of each pure fluid "
            r"(mean over all valid frames of that gas) with the number of valid frames and campaigns, and "
            r"the expanded uncertainty ($k = 2$, twice the standard deviation of the per-frame purity).}",
           r"  \begin{tabular}{lccc}", r"    \toprule",
            r"    gas & purity & $N_{\rm frames}$ ($N_{\rm campaigns}$) & $U$ ($k = 2$) \\", r"    \midrule"]
    for r in gases.itertuples():
        tex.append(f"    {r.gas} & {r.purity:.4f} & {n(r.frames)} ({r.days}) & {r.U_k2:.3f} \\\\")
    tex += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}"]
    with open(os.path.join(OUT_DIR, 'si_purity_tables.tex'), 'w') as f:
        f.write("\n".join(tex) + "\n")
    print(days.round(4).to_string(index=False))
    print(gases.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
