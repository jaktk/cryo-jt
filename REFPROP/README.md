# `REFPROP` directory

This directory contains the only non-standard part of REFPROP 10: the
mixture models in [`FLUIDS/HMX.BNC`](FLUIDS/HMX.BNC), with the reducing
parameters and departure functions for the binary mixtures:

- helium-neon (model `BHN`), Tkaczuk *et al.* (2020);
- helium-argon (model `BEA`), Tkaczuk *et al.* (2020);
- neon-argon (model `BNA`), Tkaczuk *et al.* (2020);
- helium-nitrogen (model `BEG`), a custom model from the thesis of
  J. Tkaczuk (unpublished), used for the He-N2 measurements.

## Merging into a local REFPROP installation

The analysis can only compute He-Ne and He-N2 properties once these models
are in the `HMX.BNC` file of the local REFPROP installation
(`FLUIDS/HMX.BNC` under the REFPROP root). Back up that file, then:

1. In its `#BNC` block, replace the entries for the pairs Helium/Neon,
   Helium/Argon, Neon/Argon, and Helium/Nitrogen with the four entries of the
   `#BNC` block of this repository's `FLUIDS/HMX.BNC` (each entry runs from
   its `?Name/Name` line to the next `!` line).
2. Copy the four `#MXM` blocks (`BHN`, `BEA`, `BNA`, `BEG`) into the same
   file, after its other `#MXM` blocks and before the `@END` line.

[`src/data_analysis/FluidProps.py`](../src/data_analysis/FluidProps.py) finds
REFPROP through the `RPPREFIX` environment variable, or a default
per-platform location when it is not set.
