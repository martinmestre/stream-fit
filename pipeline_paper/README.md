# Pipeline of the paper about the GD-1 stream on a MW with a fermionic DM halo

The results in the paper are a consequence of analyzing the ouputs of following pipeline.
All the plots are placed inside the directory "paper_plots".

## Fit the initial conditions (IC) of the orbit for the NFW-MW model
```
python fit_data_I-M-GaiaDR2_to_MWPot2014wGalpy.py
```

Output: "observable_orbit_NFW-MW.txt", "param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt"

## Fit the parameters of the Fermionic-MW model using the IC from NFW-MW
```
python fit_pot_from_IbataPolysGaiaDR2-data_chi2full.py
```

First run (can be skipped) edit like this:
```
bounds = ((35, 40), (25, 30), (1.1e-5, 1.4e-5))
opt = optimize.differential_evolution(chi2_full, bounds, args=(ener_f, ic, r_sun),
                                      strategy='best2bin', maxiter=200, popsize=200, tol=5.0e-8,
                                      atol=0.0, disp=True, polish=True, workers=-1)
```
The second run is to improve the solution; edit like this:
```
bounds = ((35, 37), (26, 28), (1.2e-5, 1.3e-5))
opt = optimize.differential_evolution(chi2_full, bounds, args=(ener_f, ic, r_sun),
                                      strategy='best2bin', maxiter=200, popsize=200, tol=5.0e-8,
                                      atol=0.0, disp=True, polish=True, workers=-1)
```

Output: "param_fit_pot_from_IbataPolysGaiaDR2-data_chi2full.txt"

## Make the plot for both MW model orbits in observable space

```
julia plot_observables.jl
```
Output: "observables_xxx.pdf"

## Compute the $\chi^2_{\rm{Stream}}$ function for fixed values of $(\epsilon, \beta_0)$ and for a grid in $(\theta_0,\omega_0)$ plane
```
python grid_chi2stream.py
```

To run in a large/fast cluster.

Output: "chi2stream_beta0_1.258e-05.txt"

## Compute the $\chi^2_{\rm{Stream}}$ function for fixed values of $(\epsilon, \beta_0)$ and for a grid in $(\theta_0,\omega_0-h(\theta_0))$ plane (zoomed in version).
```
python grid_chi2stream_tilted.py
```

To run in a large/fast cluster.

Output: "chi2stream_tilted_beta0_1.258e-05.txt"

## Make $\chi^2_{\rm{Stream}}$ function plot

```
julia plot_chi2stream.jl
```

Output: "chi2stream_contourf.pdf", "chi2stream_tilted_contourf.pdf"

## Compute and plot rotation curves: observations, Fermionic-MW and NFW-MW.

```
julia rotation_curves.jl
```

Output: "rotation_curves.pdf"

## Compute Fermionic-MW solutions sequencially for $\epsilon\in [56,370]$ keV

```
julia optim_polish_chi2full.jl
```

Output: "param_optim_polish_chi2full.txt"
