"""Module to share variables."""
import numpy as np

# Numerical constants
deg2rad = np.pi/180.0

# Splines
phi_1_min = 0
phi_1_max = 0
phi_2_spl = 0
d_hel_spl = 0
v_hel_spl = 0
mu_1_spl = 0
mu_2_spl = 0
mu_ra_spl = 0
mu_dec_spl = 0
ra_min = 0
ra_max = 0
dec_spl = 0
r_aux = 0

# Model parameters
g = 4.300923924e-6

# NFW (halo)
H_0 = 0.0674
rho_c_nfw = 3.0*H_0**2/(8.0*np.pi*g)
c = 12.
m_200 = 0.97e12
r_200 = (3.0*m_200/(800.0*np.pi*rho_c_nfw))**(1.0/3.0)
a_nfw = r_200/c

# Logarithmic (halo)
v_c = 220.0
q_phi = 0.9

# Plummer sphere (bulge)
M_plu = 0.5e10
r_plu = 0.23

# Miyamoto-Nagai (disk)
M_mn = 6.8e10  # Msun
a_mn = 3.0     # kpc
b_mn = 0.28    # kpc

# Exponential spheroid (inner bulge)
a_exps_in = 0.0038
rho_exps_in = 3.6e13

# Exponential spheroid (main bulge)
a_exps_ma = 0.12
rho_exps_ma = 1.9e11

# Exponential flat disk (disk)
a_efd = 3.0*0.93
m_efd = 4.4e10*1.2
sigma_efd = m_efd/(2.0*np.pi*a_efd*a_efd)
