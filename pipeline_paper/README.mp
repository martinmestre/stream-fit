# Pipeline of the paper about the GD-1 stream on a MW with a fermionic DM halo.

In order to reproduce the paper results run the following pipeline:

python fit_data_I-M-GaiaDR2_to_MWPot2014wGalpy.py
output: "observable_orbit_NFW-MW.txt", "param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt"

python fit_pot-slice_from_IbataPolysGaiaDR2-data.py
output: "param_fit_pot-slice_from_IbataPolysGaiaDR2-data.txt"

julia optim_GR.jl
output: