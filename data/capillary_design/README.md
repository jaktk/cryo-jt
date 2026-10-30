# Capillary design data

Precomputed mass-flow data from the capillary Darcy-Weisbach calculations used
when sizing the expansion capillary. It documents the longer 10 m capillary
that was considered and abandoned; no script in this repository reads or
regenerates it.

File-name convention: `D{D_int_mm}_L{L_m}_x{x_He}_p{p_in_MPa}_T{T_in_K}.csv`.
Columns:

| Column | Unit | Description |
| --- | --- | --- |
| `mf/g_per_s` | g s⁻¹ | mass-flow rate through the capillary |
| `p2/MPa` | MPa | outlet (downstream) pressure |
| `T2/K` | K | isenthalpic outlet temperature, T₂ₕ |

Currently archived: `D0.4_L10_x0.501_p13_T80.csv`, equimolar ⁴He-Ne, 80 K
inlet, 13 MPa inlet, 0.4 mm × 10 m capillary.
