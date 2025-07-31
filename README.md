# Measurements of the Joule-Thomson coefficent in cryogenic fluids

This repository contains raw data, processed data, and code used and obtain from measurements of the Joule-Thomson coefficient in pure gases and mixtures at temperature between 65 K and 100 K.

The measurements from this repository validate the earlier-developped equations of state in the single-phase region. Indirect measurements are first acquired for pure fluids, allowing for the experiment validation and then for mixtures, providing new results to the study. The expanded relative standard uncertainty is calculated along with the Monte Carlo analysis for the combined uncertainty. The impact of the composition uncertainty on the Joule-Thomson coefficient is quantified for mixtures using the Monte Carlo simulations.

## Experimental setup

All the nomenclature in this repository, the raw and processed data files in particular, follow the naming convention from the P&ID digaram below.

![pfd](img/PID.png "Piping and Instrumentation Diagram of the experimental setup for the Joule-Thomson coefficient")

The most imporant values for analysing the Joule-Thomson coefficient are:
- temperature from TT101 (upstream) and TT102 (downstream) thermometers;
- pressure from PT101 (upstream) and PT102 (downstream) thermometers.

They measure the inlet and outlet parameters for the isenthalpic expansion.

## Raw data

`raw-data` contains raw measurement data stored in CSV files. It contains:
- the pressure-temperature pairs used for indirect measurements of the Joule-Thomson coefficent in pure fluids (nitrogen, argon, helium-4);
- the pressure-temperature-composition values used for indirect measurements of the Joule-Thomson coefficient in fluid mixtures (helium-neon, helium-nitrogen).

| Variable    | Unit       | Description                                                                  |
| ----------- | ---------- | ---------------------------------------------------------------------------- |
| Date        | yyyy-mm-dd | Measurement date                                                             |
| Time        | hh:mm:ss   | Measurement time                                                             |
| PT102       | MPa        | Pressure measured after the isenthalpic expansion                            |
| PT101       | MPa        | Pressure measured before the isenthalpic expansion                           |
| TT009       | K          | Themal shield highest temperature (point of the highest radiation heat loss) |
| TT010       | K          | Temperature at the outlet from the cryostat (before heater Q011)             |
| TT008       | K          | Cold head temperature                                                        |
| TT101       | K          | Temperature measured before the isenthalpic expansion                        |
| RT101       | $\Omega$   | Resistance of the TT101 Cernox temperature sensor                            |
| TT102       | K          | Temperature measured after the isenthalpic expansion                         |
| RT102       | $\Omega$   | Resistance of the TT102 Cernox temperature sensor                            |
| Q008        | W          | Power of the cold head heater (automaticaly regulated)                       |
| x1          | -          | 1st gas concentration (for mixture measurements)                             |
| x2          | -          | 2nd gas concentration (for mixture measurements)                             |
| impurity    | -          | Gas impurity defined as $1 - purity$ (for single component measurements)     |
| err(x)      | -          | Concentration/purity measurement uncertainty                                 |

The measurements for two dates do not contain values for Q008: 2020-11-27 and 2020-12-14.

The measurements are grouped in files by collection date. The Joule-Thomson coefficient was measured for the fluids below:

| Date       | Fluid            | Average purity | Average molar composition |
| ---------- | ---------------- | -------------- | ------------------------- |
| 2020-11-27 | ${\rm N_2}$      | 0.9987         | -                         |
| 2020-12-09 | ${\rm N_2 - He}$ | -              | 0.8497/0.1503             |
| 2020-12-10 | ${\rm N_2 - He}$ | -              | 0.8909/0.1091             |
| 2020-12-11 | ${\rm N_2 - He}$ | -              | 0.5641/0.4359             |
| 2020-12-14 | ${\rm N_2}$      | 0.9904         | -                         |
| 2020-12-15 | ${\rm Ar}$       | 0.9990         | -                         |
| 2020-12-16 | ${\rm N_2}$      | 0.9977         | -                         |
| 2020-12-17 | ${\rm N_2}$      | 0.9986         | -                         |
| 2021-01-12 | ${\rm He}$       | 0.9967         | -                         |
| 2021-01-13 | ${\rm He}$       | 0.9904         | -                         |
| 2021-01-14 | ${\rm He - Ne}$  | -              | 0.7148/0.2852             |
| 2021-01-15 | ${\rm He - Ne}$  | -              | 0.7385/0.2615             |
| 2021-01-25 | ${\rm N_2}$      | 0.9977         | -                         |
| 2021-01-26 | ${\rm He}$       | 0.9993         | -                         |
| 2021-01-27 | ${\rm He - Ne}$  | -              | 0.3464/0.6536             |
| 2021-01-28 | ${\rm He - Ne}$  | -              | 0.4084/0.5916             |
| 2021-03-03 | ${\rm He - Ne}$  | -              | 0.2219/0.7781             |
| 2021-03-04 | ${\rm N_2}$      | 0.9987         | -                         |

## Code

The `src` directory consists of the following sub-directories:
- `data-acquisition` - code running data acquisition and real-time visualization of collected data
- `data-analysis` - code for data wrangling and analysis of results

	