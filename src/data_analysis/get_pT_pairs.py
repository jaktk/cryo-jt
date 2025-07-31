import os
import numpy as np
import pandas as pd
import json
from get_git_root import get_git_root


git_root = get_git_root(os.getcwd())

with open(os.path.join(git_root, "data", "time_stamps.json")) as f:
    data = json.load(f)

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
                # calculate mean, standard deviation, and standard expanded uncertainty (k=1.96)
                # for all columns but date and time
                mean_mp[col] = [mp[col].mean()]
                mean_mp[f"{col}_STD"] = [mp[col].std(skipna=True)]
                mean_mp[f"{col}_EXP_UNC"] = [2 * mean_mp[f"{col}_STD"][0] / len(mp[col])]
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
    
    # save selected measurements for a single isenthalp to csv
    measurement_isenthalp.to_csv(os.path.join(git_root,
                                              "data",
                                              "derived_data",
                                              fname),
                                 index=False)
    print("Done.")
