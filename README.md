# stream-fit
Fitting orbits, streams and Milky Way potentials to model a large range of observables including stellar streams, rotation curve and the stars at the Galactic centre.

## Pipeline of the paper about the GD-1 stream on a MW with a fermionic DM halo.

In order to reproduce the paper results go to directory "/pipeline_paper/" and run the following pipeline:

python fit_data_I-M-GaiaDR2_to_MWPot2014wGalpy.py
output: "observable_orbit_NFW-MW.txt", "param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt"

python fit_pot-slice_from_IbataPolysGaiaDR2-data.py
output: "param_fit_pot-slice_from_IbataPolysGaiaDR2-data.txt"

julia optim_GR.jl
output:
