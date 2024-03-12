"""
Grid of potentials and their likelihoods.

Author: Martín Mestre.
"""


import numpy as np
import itertools as iter
import matplotlib.pyplot as plt
import multiprocessing as mp
import h5py
import stream


def worker(theta_0, d_theta, beta_0):
    """Parallel task."""
    print(theta_0, d_theta, beta_0)
    return stream.chi2_stream(theta_0, d_theta, beta_0, ener_f, ic, r_sun)



# Parameters
n_beta = 1
n_grid = 500
r_sun = 8.122  # kpc   (Gravity Collaboration 2018.)
ener_f = 56.0  # keV
bounds = ((34.0, 38.0), (25, 30))  # For (theta_0, W_0-theta_0)
ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
ic = np.loadtxt(ic_file)
n_job = mp.cpu_count()
print('n_cpu = ', n_job)
param_file = "serafin/sol_optim_pot_m{:2d}.txt".format(int(ener_f))
_a, _b, beta_b = np.loadtxt(param_file)
beta_lim = np.linspace(beta_b, beta_b+1, n_beta)


if __name__ == "__main__":

    for beta in beta_lim:
        chi2stream_file = 'dirac/chi2stream_beta0_{:.3e}.txt'.format(beta)

        args_iter = list(iter.product(np.linspace(bounds[0][0], bounds[0][1], n_grid),
                                      np.linspace(bounds[1][0], bounds[1][1], n_grid),
                                      beta*np.ones(1)))

        with mp.Pool(processes=n_job) as pool:
            # starmap is supposed to give the result ordered.
            results = pool.starmap(worker, iterable=args_iter)
            chi = list(results)
        print('results=', chi)


        # Open likelhood file
        f = open(chi2stream_file, 'wt')

        for i in range(0, len(args_iter)):
            args = args_iter[i]
            print('args =', args)
            theta_0, d_theta, beta_0 = args
            f.write("{}  {}  {}  {}".format(theta_0, d_theta, beta_0, chi[i]))
            f.write("\n")

            i += 1
        f.close()
