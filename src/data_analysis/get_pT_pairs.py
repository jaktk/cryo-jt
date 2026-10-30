import os
import numpy as np
import pandas as pd
import json
from repo_root import REPO_ROOT


git_root = REPO_ROOT
DRIFT_INTERVAL = 60  # s before each window over which the drift of the outlet state is fitted (measurement_uncertainty.DRIFT_INTERVAL)
DRIFT_COLUMNS = ["PT102/MPa", "TT102/K"]

with open(os.path.join(git_root, "data", "raw_data", "extraction_windows.json")) as f:
    data = json.load(f)

written = {}
for isenthalp_stamps in data["isenthalps"]:
    print(f"Extracting points from isenthalp {isenthalp_stamps['id']}...", end=" ")
    _stamp = ""
    measurement_isenthalp = pd.DataFrame()
    for stamp in isenthalp_stamps["stamps"]:
        if _stamp != stamp:
            # read correct dataframe only if it is not in current memory
            df = pd.read_csv(os.path.join(git_root,
                                          "data",
                                          "raw_data",
                                          f"{stamp['date/yyyy-mm-dd']}.csv"),
                             sep = ",",
                             engine = "python")
            _stamp = stamp
        # find measurement points (mp) in csv files containing time stamps specified in json file
        mp = df.loc[(df['time/hh:mm:ss'] >= stamp["time1/hh:mm:ss"])
                    & (df['time/hh:mm:ss'] <= stamp["time2/hh:mm:ss"])]
        mean_mp = {}
        for col in mp.columns:
            if col not in ["date/yyyy-mm-dd", "time/hh:mm:ss"]:
                # calculate mean, standard deviation, and expanded uncertainty of the
                # mean (k=2) for all columns but date and time. The expanded
                # uncertainty uses the standard error of the mean, std / sqrt(N).
                mean_mp[col] = [mp[col].mean()]
                mean_mp[f"{col}_STD"] = [mp[col].std(skipna=True)]
                mean_mp[f"{col}_EXP_UNC"] = [2.0 * mean_mp[f"{col}_STD"][0] / np.sqrt(len(mp[col]))]
        # drift of the outlet state: slope (per s) of a straight line fitted
        # over the DRIFT_INTERVAL seconds before the window
        seconds = pd.to_timedelta(df['time/hh:mm:ss']).dt.total_seconds()
        start = pd.to_timedelta(stamp["time1/hh:mm:ss"]).total_seconds()
        before = df.loc[(seconds >= start - DRIFT_INTERVAL) & (seconds < start)]
        t_before = seconds.loc[before.index]
        for col in DRIFT_COLUMNS:
            mean_mp[f"{col}_DRIFT"] = [np.polyfit(t_before, before[col], 1)[0]]
        measurement_isenthalp = pd.concat([measurement_isenthalp,
                                           pd.DataFrame.from_dict(mean_mp)])
    measurement_isenthalp.reset_index(drop=True, inplace=True)
    
    # set unique name to new csv file with selected measurement data
    Tin = round(measurement_isenthalp["TT101/K"].mean())
    pin = round(measurement_isenthalp["PT101/MPa"].mean())
    if type(isenthalp_stamps['fluid']) is list:
        fluid_name = "-".join(isenthalp_stamps['fluid'])
    else:
        fluid_name = isenthalp_stamps['fluid']
    fname = f"{fluid_name}_{Tin}K_{pin}MPa.csv"
    if fname in written:
        raise ValueError(f"isenthalps {written[fname]} and {isenthalp_stamps['id']} "
                         f"both map to {fname}; exclude one in extraction_windows.json")
    written[fname] = isenthalp_stamps['id']
    
    # save selected measurements for a single isenthalp to csv
    measurement_isenthalp.to_csv(os.path.join(git_root,
                                              "data",
                                              "derived_data",
                                              "p_T_pairs",
                                              fname),
                                 index=False)
    print("Done.")
