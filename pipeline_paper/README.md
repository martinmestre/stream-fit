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
bounds = ((35, 40), (25, 30), (1.0e-5, 1.5e-5))
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

## Improve the ICs for the Fermionic-MW model
```
python fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.py
```

Output: "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"


## Re-fit the parameters for the Fermionic-MW model using NOMAD solver in small box,using Distributed.jl parallel scheme in a SLURM cluster environment.

```
sbatch -N 3 --ntasks-per-node=32 --partition=batch -o optim_pot_m56.out optim_pot.jl 1
```
having set before:

```
const lb_g = [[35.8, 27.0, 1.2e-5], [36., 27., 1.2e-5], [37., 28., 5.0e-5],
              [38., 29., 3.5e-4], [40., 29., 1.3e-3], [43., 29.6, 3.0e-3]]
const ub_g = [[36.3, 27.6, 1.3e-5], [40., 31., 1.0e-4], [41., 32., 1.0e-3],
              [42., 32., 3.0e-3], [44., 32., 4.0e-3], [47., 36., 1.0e-2]]

const n_grid = 4
```


Output: "sol_optim_pot_m56.txt"  "chi2_optim_pot_m56.txt"


## Make the plot for both MW model orbits in observable space

```
julia plot_observables.jl
```
Output: "observables_xxx.pdf"

## Compute the $\chi^2_{\rm{Stream}}$ function for fixed values of $(\epsilon=56, \beta_0\approx1.254\times10^{-5})$ and for a grid in $(\theta_0,\omega_0)$ plane
```
python grid_chi2stream.py
```

To run in a large/fast cluster.

Output: "chi2stream_beta0_1.254e-05.txt"

## Compute the $\chi^2_{\rm{Stream}}$ function for fixed values of $(\epsilon, \beta_0)$ and for a grid in $(\theta_0,\omega_0-h(\theta_0))$ plane (zoomed in version).
```
python grid_chi2stream_tilted.py
```

To run in a large/fast cluster.

Output: "chi2stream_tilted_beta0_1.254e-05.txt"

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

## Compute Fermionic-MW solutions for any $\epsilon$ using Distributed.jl parallel scheme in a SLURM cluster environment.

In the example below replace the $M \in [1:5]$ according to
the index of the fermion mass $m=[56, 100, 200, 300, 360]$.
```
sbatch -N 3 --ntasks-per-node=64 --partition=multi -o optimo_pot_imM.out optim_pot.jl M
```

Output: "sol_optim_pot_mF.txt" where $F\in m$.



## Plot the density profiles for $\epsilon=56, 100, 200, 300, 360$ keV$

```
julia plot_density_profiles.jl
```

Output: "density_profiles.pdf"