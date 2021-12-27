"""Test example for RAR model."""

import numpy as np
import model_def_units as model_def
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
import astropy.units as u
import time

# Model evaluation
ener_f = 56.0*u.keV
# theta_0, W_0, beta_0 = 3.77780827e+01, 6.63468885e+01, 1.20446329e-05
theta_0 = 3.623511473208260014e+01
W_0 = theta_0 + 2.745737545624503895e+01
beta_0 = 1.1977e-5
param = [ener_f, theta_0, W_0, beta_0]
start_time = time.time()
r, mass, nu = model_def.model(param)
execution_time = (time.time()-start_time)
print('time=', execution_time)

# Density profile
mass_spline = InterpolatedUnivariateSpline(r.value, mass.value, k=4)  # Allows easy computation of derivatives


def rho_spline(r):
    """Density function, afterwards."""
    deriv = mass_spline.derivative(1)
    return deriv(r)/(4.0*np.pi*r*r)


r = r.value
mass = mass.value
rho = rho_spline(r)

np.savetxt('test.out', np.column_stack((r, mass, nu)))

# Plots
fig = plt.figure(figsize=(10, 7))
font = {"size": 15}
plt.rc('font', **font)
plt.scatter((4.2538e-7, 12.e0, 40.e0), (3.5e6, 3.625708546037865e+10, 2.273879460291010e+11))
plt.plot(r, mass)
# plt.xscale('log')
# plt.yscale('log')
# plt.ylim(0, 1.e12)
plt.xlabel(r'$r$ [kpc]')
plt.ylabel(r'$M$ [M$_\odot$]')
plt.xlim(1.e-10, r[-1])
plt.show()

fig = plt.figure(figsize=(10, 7))
font = {"size": 15}
plt.rc('font', **font)
plt.plot(r, rho)
plt.xscale('log')
plt.yscale('log')
# plt.ylim(1.e+4,1.e35)
# plt.xlim(1.e-5,1.e-3)
plt.xlabel(r'$r$ [kpc]')
plt.ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$]')
plt.show()

fig = plt.figure(figsize=(10, 7))
font = {"size": 15}
plt.rc('font', **font)
plt.plot(r, nu)
plt.xscale('log')
# plt.yscale('log')
# plt.ylim(0,1.e12)
plt.xlabel(r'$r$ [kpc]')
plt.ylabel(r'$\nu$')
plt.xlim(1.e-10, 1.e2)

plt.show()
