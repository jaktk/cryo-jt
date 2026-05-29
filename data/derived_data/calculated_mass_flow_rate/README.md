# `data/derived_data/calculated_mass_flow_rate/`

Capillary-sizing design calculations: downstream pressure as a function of
mass flow rate for the chosen 6 m, 0.365 mm ID stainless-steel capillary, for
each inlet condition used to plan the campaigns. These are calculations, not
measurements; the calculation is described in §2.3 of the paper.

File names encode the capillary geometry, mixture composition, and inlet
condition:

`Helium-Neon_D0.365_L6_x0.501_p10_T80.csv`

| Token | Meaning |
| --- | --- |
| `Fluid` | Fluid or mixture (`Helium-Neon`, `Helium-Nitrogen`, etc.) |
| `D<value>` | Capillary internal diameter, mm |
| `L<value>` | Capillary length, m |
| `x<value>` | Mole fraction of the first mixture component |
| `p<value>` | Inlet pressure, MPa |
| `T<value>` | Inlet temperature, K |

Column definitions are in
[`../../metadata/codebook.csv`](../../metadata/codebook.csv).
