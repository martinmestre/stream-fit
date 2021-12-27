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
    log_like = stream.log_likelihood(theta_0, d_theta, beta_0, ener_f, ic, r_sun)
    return log_like


def save_hdf5(result, file):
    """Save with hdf5. Starmap gives the result ordered."""
    f = h5py.File(file, 'w')

    for idx, data in enumerate(result):
        group_name = 'models_{:02d}'.format(idx)
        dset = f.create_group(group_name)
        dset.attrs['theta_0'] = result[idx][0]
        dset.attrs['d_theta'] = result[idx][1]
        dset.attrs['beta_0'] = result[idx][2]
        dset.attrs['loglike'] = result[idx][3]
        len_r = len(result[idx][4])
        r = dset.create_dataset('r', (len_r,), dtype='f')
        mass = dset.create_dataset('mass', (len_r,), dtype='f')
        nu = dset.create_dataset('nu', (len_r,), dtype='f')
        r[:] = result[idx][4]
        mass[:] = result[idx][5]
        nu[:] = result[idx][6]
        print('r = ', r)

    print("Closing HDF5 file")

    f.close()

    return None


# Parameters
n_beta = 20
n_grid = 500
r_sun = 8.122  # kpc   (Gravity Collaboration 2018.)
ener_f = 56.0  # keV
bounds = ((35, 37), (26, 28), (1.0e-5, 1.3e-5))  # For (theta_0, W_0-theta_0, beta_0)
ic = np.array([1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
              -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01])
n_job = mp.cpu_count()
print('n_cpu = ', n_job)
beta_lim = np.linspace(1.0e-5, 1.5e-5, n_beta)


if __name__ == "__main__":

    for beta in beta_lim:
        likelihood_file = 'likelihood_beta0_{:.2e}.txt'.format(beta)
        h5_file = 'model_beta0_{:.2e}.hdf5'.format(beta)

        args_iter = list(iter.product(np.linspace(bounds[0][0], bounds[0][1], n_grid),
                                      np.linspace(bounds[1][0], bounds[1][1], n_grid),
                                      beta*np.ones(1)))

        with mp.Pool(processes=n_job) as pool:
            # starmap is supposed to give the result ordered.
            results = pool.starmap(worker, iterable=args_iter)
            ll = list(results)
        print('results=', ll)

        # Save models to a HDF5 file
        save_hdf5(results, h5_file)

        # Open likelhood file
        f = open(likelihood_file, 'wt')

        for i in range(0, len(args_iter)):
            args = args_iter[i]
            print('args =', args)
            theta_0, d_theta, beta_0 = args
            f.write("{}  {}  {}  {}".format(theta_0, d_theta, beta_0, ll[i][3]))
            f.write("\n")

            i += 1
        f.close()
