# Pipeline of the paper about the GD-1 stream on a MW with a fermionic DM halo.

The results in the paper are a consequence of analyzing the ouputs of following pipeline:

### python fit_data_I-M-GaiaDR2_to_MWPot2014wGalpy.py

output: "observable_orbit_NFW-MW.txt", "param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt"

### python fit_pot_from_IbataPolysGaiaDR2-data_chi2full.py

First run (can be skipped) edit like this:
```
bounds = ((35, 40), (25, 30), (1.1e-5, 1.4e-5))
opt = optimize.differential_evolution(chi2_full, bounds, args=(ener_f, ic, r_sun),
                                      strategy='best2bin', maxiter=200, popsize=200, tol=5.0e-8,
                                      atol=0.0, disp=True, polish=True, workers=-1)
```
The second run is to polish the solution; edit like this:
```
bounds = ((35, 37), (26, 28), (1.2e-5, 1.3e-5))
opt = optimize.differential_evolution(chi2_full, bounds, args=(ener_f, ic, r_sun),
                                      strategy='best2bin', maxiter=200, popsize=200, tol=5.0e-8,
                                      atol=0.0, disp=True, polish=True, workers=-1)
```
output: "param_fit_pot_from_IbataPolysGaiaDR2-data_chi2full.txt"


### julia optim_polish_chi2full.jl

output: "param_optim_polish_chi2full.txt"
