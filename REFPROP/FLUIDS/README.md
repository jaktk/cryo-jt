# `REFPROP/FLUIDS/`

Custom mixture-model coefficients for REFPROP: the He-Ne, He-Ar, and Ne-Ar
binary mixtures developed in the prior work of
[Tkaczuk *et al.* 2020](https://doi.org/10.1063/1.5142275), and a custom
He-N2 model (J. Tkaczuk, thesis, unpublished).
The file `HMX.BNC` here must be merged into REFPROP's own `FLUIDS/HMX.BNC`
before the analysis pipeline can compute mixture properties; the steps are in
the parent [`../README.md`](../README.md#merging-into-a-local-refprop-installation).
