"""
Fit GD-1 stream data varying a slice of the RAR's parameter but keeping fixed the orbit initial conditions.

Computation of the orbit model for RAR + barionic component.
Using Ibata+20 polinomial data-fits.
Optimization to fit the best RAR parameters.

Author: Martín Mestre.
"""

import numpy as np
import matplotlib.pyplot as plt
import config as cfg
import astropy.coordinates as coord
import GD1Koposov10_class as GD1_class
import pandas as pd
import potential_classes as pot
from astropy import units as u
from scipy.integrate import solve_ivp
from scipy import optimize
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from scipy.signal import argrelextrema

def accel_mw(pot_list, x, y, z):
    """Acceleration function generated by the MW total potential."""
    bulge, thin, thick, halo = pot_list
    return bulge.accel(x, y, z)+thin.accel(x, y, z)+thick.accel(x, y, z)+halo.accel(x, y, z)


def symp_grad_mw(t, w, pot_list):
    """Symplectic gradient generated by MW total potential."""
    x = w[0]
    y = w[1]
    z = w[2]
    px = w[3]
    py = w[4]
    pz = w[5]
    acceleration = accel_mw(pot_list, x, y, z)
    return [px, py, pz, acceleration[0], acceleration[1], acceleration[2]]


def grad_mw(pot_list, x, y, z):
    """Minus the acceleration."""
    return -accel_mw(pot_list, x, y, z)


def rot_vel_mw(pot_list, r):
    """Circular velocity."""
    return np.sqrt(r*(grad_mw(pot_list, r, 0, 0)[0]))

def v_circ_GR(r, halo):
    """GR circula velocity for the rar model. Better to use this one instead of the Newtonian."""
    c = 2.997925e5  # km s^-1
    return np.sqrt(0.5*c*c*r*halo.dnu_wrap(r))


def get_core_GR(halo):
    """Compute the core mass and radius."""
    arg_max = argrelextrema(v_circ_GR(halo.r_s, halo), np.greater)  # The approximation is Newtonian
    arg_first_max = arg_max[0][0]
    r_cand = halo.r_s[arg_first_max]
    bounds = np.array([r_cand*0.5, r_cand*1.5])
    r_core = optimize.fminbound(lambda r: -v_circ_GR(r, halo),
                                bounds[0], bounds[1], xtol=0.5e-12, maxfun=1000)
    m_core = halo.mass_wrap(r_core)
    print('r_cand = ', r_cand, ' ---- r_core = ', r_core)
    return r_core, m_core

def pot_model(ener_f, theta_0, W_0, beta_0):
    """Potential model including barions and dark matter halo."""
    # Set barionic potentials
    M_gal = 2.32e7  # Msun
    bulge = pot.Plummer(460.0*M_gal, 0.3)
    thin = pot.MiyamotoNagai(1700.0*M_gal, 5.3, 0.25)
    thick = pot.MiyamotoNagai(1700.0*M_gal, 2.6, 0.8)
    # Set halo potential
    param = np.array([ener_f, theta_0, W_0, beta_0])
    halo = pot.RAR(param)
    return bulge, thin, thick, halo


def orbit_model(alpha, delta, distance, mu_alpha, mu_delta, v_los, pot_list, r_sun):
    """
    Orbit definition.

    The orbit_model here defined works with sky coordinates
    at input and sky-cartesian at output.
    """
    bulge, thin, thick, halo = pot_list

    # Autoconsistent velocity of LSR
    v_circ_sun = rot_vel_mw(pot_list, r_sun)

    # Transformation to galactocentric coordinates
    sky_coord = coord.ICRS(ra=alpha*u.degree, dec=delta*u.degree,
                           distance=distance*u.kpc,
                           pm_ra_cosdec=mu_alpha*np.cos(delta*u.degree)*u.mas/u.yr,
                           pm_dec=mu_delta*u.mas/u.yr,
                           radial_velocity=v_los*u.km/u.s)
    galcen_distance = r_sun*u.kpc
    v_sun = coord.CartesianDifferential([11.1, v_circ_sun+12.24, 7.25]*u.km/u.s)
    z_sun = 0.0*u.kpc
    frame = coord.Galactocentric(galcen_distance=galcen_distance, galcen_v_sun=v_sun, z_sun=z_sun)
    galac_coord = sky_coord.transform_to(frame)

    w_0 = np.zeros(6)
    w_0[:3] = [galac_coord.x/u.kpc, galac_coord.y/u.kpc, galac_coord.z/u.kpc]
    w_0[3:] = [galac_coord.v_x/(u.km/u.s), galac_coord.v_y/(u.km/u.s), galac_coord.v_z/(u.km/u.s)]

    # ODE integration
    unit_t = 0.977792221680356   # Gyr
    time_span_s2 = 0.2/unit_t
    t_0 = 0.0/unit_t
    n_steps = 1000
    t_back = np.linspace(t_0, -time_span_s2, n_steps+1)
    t_forw = np.linspace(t_0, time_span_s2, n_steps+1)
    sol_back = solve_ivp(symp_grad_mw, [t_0, -time_span_s2], w_0, t_eval=t_back, args=[pot_list],
                         method='DOP853', rtol=5.0e-14, atol=0.5e-14)
    sol_forw = solve_ivp(symp_grad_mw, [t_0, time_span_s2], w_0, t_eval=t_forw, args=[pot_list],
                         method='DOP853', rtol=5.0e-14, atol=0.5e-14)
    # t = np.concatenate([sol_back.t,sol_forw.t])  Not used
    y = np.concatenate([sol_back.y, sol_forw.y],  axis=1)
    y = np.delete(y, 0, axis=1)  # Remove duplicated column

    # Transformation to GD-1 frame of coordinates (\phi_1, \phi_2)
    galac_coord = coord.Galactocentric(x=y[0]*u.kpc, y=y[1]*u.kpc, z=y[2]*u.kpc,
                                       v_x=y[3]*u.km/u.s, v_y=y[4]*u.km/u.s, v_z=y[5]*u.km/u.s,
                                       galcen_distance=galcen_distance, galcen_v_sun=v_sun, z_sun=z_sun)
    gd1_coord = galac_coord.transform_to(GD1_class.GD1Koposov10)
    phi_1 = gd1_coord.phi1
    phi_2 = gd1_coord.phi2
    d_hel = gd1_coord.distance
    v_hel = gd1_coord.radial_velocity
    # Unused block of code
    # mu_phi_1 = gd1_coord.pm_phi1_cosphi2/np.cos(phi_2)  #not used by Ibata
    # mu_phi_2 = gd1_coord.pm_phi2
    # return phi_1, phi_2, d_hel, v_hel, mu_phi_1, mu_phi_2
    # Transformation to ICRS coordinates
    icrs_coord = galac_coord.transform_to(coord.ICRS)
    mu_ra = icrs_coord.pm_ra_cosdec  # Ibata's mu_ra = pm_ra_cosdec
    mu_dec = icrs_coord.pm_dec
    return phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel, y[0], y[1], y[2], v_circ_sun


class IbaPoly:
    """
    Ibata polynomials.

    [x] = radians
    [S] = radians
    [D] = kpc
    [V] = km/s
    [MU] = mas/year
    """

    def S(self, x):
        """Sky position."""
        return 0.008367*x**3-0.05332*x**2-0.07739*x-0.02007

    def D(self, x):
        """Photometric distance."""
        return -4.302*x**5-11.54*x**4-7.161*x**3+5.985*x**2+8.595*x+10.36

    def V(self, x):
        """Radial velocity."""
        return 90.68*x**3+204.5*x**2-254.2*x-261.5

    def MU_RA(self, x):
        """Proper motion in RA (without cos(Dec)."""
        return 3.794*x**3+9.467*x**2+1.615*x-7.844

    def MU_DEC(self, x):
        """Proper motion in Declination."""
        return -1.225*x**3+8.313*x**2+18.68*x-3.95

    limit = [-90, 10]


# Observations (Ibata polynomials evaluated in a grid)
pol = IbaPoly()
Iba_sky = pd.DataFrame()
Iba_sky['phi_1'] = np.linspace(IbaPoly.limit[0], IbaPoly.limit[1], 100)
Iba_sky['phi_2'] = pol.S(cfg.deg2rad*Iba_sky['phi_1'])/cfg.deg2rad
Iba_sky['d_hel'] = pol.D(cfg.deg2rad*Iba_sky['phi_1'])
Iba_sky['v_hel'] = pol.V(cfg.deg2rad*Iba_sky['phi_1'])
Iba_sky['mu_ra'] = pol.MU_RA(cfg.deg2rad*Iba_sky['phi_1'])
Iba_sky['mu_dec'] = pol.MU_DEC(cfg.deg2rad*Iba_sky['phi_1'])

# Transform Ibata data to Galactocentric coordinates in order to know
# the galactocentricdistances.
Iba_coord_GD1 = GD1_class.GD1Koposov10(phi1=Iba_sky['phi_1']*u.degree, phi2=Iba_sky['phi_2']*u.degree)
Iba_coord_ICRS = Iba_coord_GD1.transform_to(coord.ICRS)
Iba_sky['ra'] = Iba_coord_ICRS.ra.value
Iba_sky['dec'] = Iba_coord_ICRS.dec.value
Iba_coord_ICRS = coord.ICRS(ra=Iba_sky['ra']*u.degree,
                            dec=Iba_sky['dec']*u.degree,
                            distance=Iba_sky['d_hel']*u.kpc,
                            pm_ra_cosdec=Iba_sky['mu_ra']*u.mas/u.yr,
                            pm_dec=Iba_sky['mu_dec']*u.mas/u.yr,
                            radial_velocity=Iba_sky['v_hel']*u.km/u.s)
z_sun = 0.0*u.kpc
frame = coord.Galactocentric(galcen_distance=8.0*u.kpc, z_sun=0.0*u.kpc)
Iba_coord_Gal = Iba_coord_ICRS.transform_to(frame)
Iba_d_gal = np.zeros(len(Iba_coord_Gal.x.value))
for i in range(0, len(Iba_coord_Gal.x.value)):
    Iba_d_gal[i] = np.sqrt(Iba_coord_Gal.x[i].value**2+Iba_coord_Gal.y[i].value**2
                           + Iba_coord_Gal.z[i].value**2)
print('Galactocentric distances of GD-1:')
print(Iba_d_gal)

np.savetxt('galcoords.txt', np.column_stack([Iba_coord_Gal.x.value, Iba_coord_Gal.y.value,
                                             Iba_coord_Gal.z.value]))

# Alternative computation of helicentric distance as a function of phi_1, by
# integrating the heliocentric velocity. The development comes from the two
# following definitions:
# 1) mu_1 = (d/dt)phi_1[mas/year]*cos(phi_2) = alpha * (d/dt)phi_1[rad/sec]*cos(phi_2)
#         with alpha = rad2sec / sec2year
# 2) v_hel(phi_1) = (d/dt)r_hel = (d/dphi_1)r_hel* (d/dt)phi_1
#
# (1) en (2) => (d/dphi_1)r_hel = v_hel(phi_1) / (d/dt)phi_1 =
#                               = v_hel(phi_1) /(mu_1/alpha/cos(phi_2)) =
#                               = v_hel(phi_1)*alpha*cos(phi_2)/mu_1 =
#                               = alpha*v_hel(phi_1)*cos(phi_2(phi_1))/mu_1(phi_1).

Iba_coord_GD1 = Iba_coord_ICRS.transform_to(GD1_class.GD1Koposov10)
print('Iba_coord_GD1 =', Iba_coord_GD1)
Iba_sky['mu_1'] = Iba_coord_GD1.pm_phi1_cosphi2
cm2kpc = 1.0/3.08567758e21
km2cm = 10.0**5
km2kpc = km2cm*cm2kpc
kpc2km = 3.08567758e16
rad2deg = 180.0/np.pi
deg2mas = 60.0*60.0*1000.0
year2sec = 365.25*24.0*60.0*60.0
sec2year = 1.0/year2sec
alpha = deg2mas/sec2year
gamma = kpc2km
integrand = (alpha/gamma)*np.cos(Iba_sky['phi_2']*u.deg)*Iba_sky['v_hel']/Iba_sky['mu_1']
ang = Iba_sky['phi_1']
r_hel = np.zeros(len(ang))
r_hel[0] = 10.4
for i in range(1, len(ang)):
    r_hel[i] = r_hel[0] + integrate.simps(integrand[0:i+1], ang[0:i+1], even='avg')
print('r_hel =', r_hel)
Iba_coord_ICRS = coord.ICRS(ra=Iba_sky['ra']*u.degree,
                            dec=Iba_sky['dec']*u.degree,
                            distance=r_hel*u.kpc,
                            pm_ra_cosdec=Iba_sky['mu_ra']*u.mas/u.yr,
                            pm_dec=Iba_sky['mu_dec']*u.mas/u.yr,
                            radial_velocity=Iba_sky['v_hel']*u.km/u.s)
z_sun = 0.0*u.kpc
frame = coord.Galactocentric(galcen_distance=8.0*u.kpc, z_sun=0.0*u.kpc)
Iba_coord_Gal = Iba_coord_ICRS.transform_to(frame)
Iba_r_gal = np.zeros(len(Iba_coord_Gal.x.value))
for i in range(0, len(Iba_coord_Gal.x.value)):
    Iba_r_gal[i] = np.sqrt(Iba_coord_Gal.x[i].value**2+Iba_coord_Gal.y[i].value**2
                           + Iba_coord_Gal.z[i].value**2)

# Core constraint
m_core_const = 3.5e6  # M_sun

sigma_array = np.array([0.5, 1.5, 2.0, 2.0, 10.0])

# Chi^2 function
def chi2(w_0, ener_f, beta_0, ic):
    """Chi^2 function."""
    import wrap
    theta_0 = w_0[0]
    d_theta = w_0[1]
    W_0 = theta_0 + d_theta
    pot_list = pot_model(ener_f, theta_0, W_0, beta_0)
    phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel, x, y, z, v_circ = orbit_model(
        ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], pot_list)
    cfg.phi_2_spl = interp1d(phi_1, phi_2, kind='cubic')
    cfg.d_hel_spl = interp1d(phi_1, d_hel, kind='cubic')
    cfg.v_hel_spl = interp1d(phi_1, v_hel,  kind='cubic')
    cfg.mu_ra_spl = interp1d(phi_1, mu_ra, kind='cubic')
    cfg.mu_dec_spl = interp1d(phi_1, mu_dec, kind='cubic')
    cfg.phi_1_min = np.amin(phi_1.value)
    cfg.phi_1_max = np.amax(phi_1.value)

    sum = np.zeros(5)

    y_mod = wrap.phi_2_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['phi_2']
    sigma2 = 0.5**2
    sum[0] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.d_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['d_hel']
    sigma2 = 1.5**2
    sum[1] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_ra_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_ra']
    sigma2 = 4.0
    sum[2] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_dec_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_dec']
    sigma2 = 4.0
    sum[3] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.v_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['v_hel']
    sigma2 = 100.0
    sum[4] = np.sum((y_dat-y_mod)**2 / sigma2)

    print('chi^2 =', np.sum(sum), '  v_circ = ', v_circ)
    return np.sum(sum)

def chi2_stream(w_0, ener_f, ic, r_sun):
    """Chi^2 stream function."""
    import wrap
    theta_0 = w_0[0]
    d_theta = w_0[1]
    beta_0  = w_0[2]
    W_0 = theta_0 + d_theta
    print("params = ",theta_0, d_theta, beta_0)
    pot_list = pot_model(ener_f, theta_0, W_0, beta_0)
    phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel, x, y, z, v_circ = orbit_model(
        ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], pot_list, r_sun)
    cfg.phi_2_spl = interp1d(phi_1, phi_2, kind='cubic')
    cfg.d_hel_spl = interp1d(phi_1, d_hel, kind='cubic')
    cfg.v_hel_spl = interp1d(phi_1, v_hel,  kind='cubic')
    cfg.mu_ra_spl = interp1d(phi_1, mu_ra, kind='cubic')
    cfg.mu_dec_spl = interp1d(phi_1, mu_dec, kind='cubic')
    cfg.phi_1_min = np.amin(phi_1.value)
    cfg.phi_1_max = np.amax(phi_1.value)

    sum = np.zeros(5)

    y_mod = wrap.phi_2_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['phi_2']
    sigma2 = sigma_array[0]**2
    sum[0] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.d_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['d_hel']
    sigma2 = sigma_array[1]**2
    sum[1] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_ra_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_ra']
    sigma2 = sigma_array[2]**2
    sum[2] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_dec_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_dec']
    sigma2 = sigma_array[3]**2
    sum[3] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.v_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['v_hel']
    sigma2 = sigma_array[4]**2
    sum[4] = np.sum((y_dat-y_mod)**2 / sigma2)

    halo = pot_list[3]

    r_core, mass_core = get_core_GR(halo)

    print(theta_0, d_theta, beta_0, ener_f, '--- chi2_stream = ', np.sum(sum))
    print("::: r_core = ", r_core, " ::: m_core = ", mass_core/1.e6, "x10⁶ M_sun")
    print("-----------")
    return np.sum(sum)


def chi2_core(w_0, ener_f, ic, r_sun):
    """Chi^2 core function."""
    import wrap
    theta_0 = w_0[0]
    d_theta = w_0[1]
    beta_0  = w_0[2]
    W_0 = theta_0 + d_theta
    param = np.array([ener_f, theta_0, W_0, beta_0])
    halo = pot.RAR(param)
    r_core, mass_core = get_core_GR(halo)
    chisq_core = (mass_core - m_core_const)**2/(0.01*m_core_const)**2
    print(theta_0, d_theta, beta_0, ener_f, '--- chi2_core = ', chisq_core)
    print("::: r_core = ", r_core, " ::: m_core = ", mass_core/1.e6, "x10⁶ M_sun")
    print("-----------")
    return chisq_core


def chi2_full(w_0, ener_f, ic, r_sun):
    """Chi^2 core function."""
    import wrap
    theta_0 = w_0[0]
    d_theta = w_0[1]
    beta_0  = w_0[2]
    W_0 = theta_0 + d_theta
    print("params = ",theta_0, d_theta, beta_0)
    pot_list = pot_model(ener_f, theta_0, W_0, beta_0)
    phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel, x, y, z, v_circ = orbit_model(
        ic[0], ic[1], ic[2], ic[3], ic[4], ic[5], pot_list, r_sun)
    cfg.phi_2_spl = interp1d(phi_1, phi_2, kind='cubic')
    cfg.d_hel_spl = interp1d(phi_1, d_hel, kind='cubic')
    cfg.v_hel_spl = interp1d(phi_1, v_hel,  kind='cubic')
    cfg.mu_ra_spl = interp1d(phi_1, mu_ra, kind='cubic')
    cfg.mu_dec_spl = interp1d(phi_1, mu_dec, kind='cubic')
    cfg.phi_1_min = np.amin(phi_1.value)
    cfg.phi_1_max = np.amax(phi_1.value)

    sum = np.zeros(6)

    y_mod = wrap.phi_2_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['phi_2']
    sigma2 = 0.5**2
    sum[0] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.d_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['d_hel']
    sigma2 = 1.5**2
    sum[1] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_ra_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_ra']
    sigma2 = 4.0
    sum[2] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_dec_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_dec']
    sigma2 = 4.0
    sum[3] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.v_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['v_hel']
    sigma2 = 100.0
    sum[4] = np.sum((y_dat-y_mod)**2 / sigma2)

    halo = pot_list[3]

    r_core, mass_core = get_core_GR(halo)
    # print('GR: ', r_core, mass_core)
    # r_core, mass_core = get_core(halo.r, halo.mass_spline)
    # print('Newton: ', r_core, mass_core)
    sum[5] = (mass_core - m_core_const)**2/(0.01*m_core_const)**2

    print(theta_0, d_theta, beta_0, ener_f,  '--- chi2_full = ', np.sum(sum))
    print('with ::: chi2_stream = ', np.sum(sum[0:5]), ' ::: chi2_core =', sum[5])
    print("and  ::: r_core = ", r_core, " ::: m_core = ", mass_core/1.e6, "x10⁶ M_sun")
    print("-----------")
    return np.sum(sum)


# Initial condition for the orbit

def invert_ic(u_0):
    """Transform from GD-1 coordinates to ICRS for the initial conditions."""
    w_0 = u_0
    gd1_coord = GD1_class.GD1Koposov10(phi1=u_0[0]*u.degree, phi2=u_0[1]*u.degree)
    icrs_coord = gd1_coord.transform_to(coord.ICRS)
    w_0[0] = icrs_coord.ra.value
    w_0[1] = icrs_coord.dec.value
    return w_0

# We take the initial condition from the good data from Ibata.
# k0 = 50
# u_0 = np.array([Iba_sky['phi_1'][k0], Iba_sky['phi_2'][k0], Iba_sky['d_hel'][k0],
#                 Iba_sky['mu_ra'][k0], Iba_sky['mu_dec'][k0], Iba_sky['v_hel'][k0]])
# ic = invert_ic(u_0)
# print('ic=', ic)


# Taking the initial conditions from the Galpy fit with fixed MW2014 potential and Ibata polinomials-data.
ic_file = "param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt"
ic = np.loadtxt(ic_file)
# Parameters
param_file = "param_fit_pot_from_IbataPolysGaiaDR2_chi2full.txt"
r_sun = 8.122  # Gravity Collaboration (2018)
ener_f = 56.0  # keV
beta_0 = 1.20e-5


#Optimization
# bounds = ((35, 37), (26, 28), (1.2e-5, 1.3e-5))
# bounds = ((35, 40), (25, 30), (1.0e-5, 1.5e-5))
# opt = optimize.differential_evolution(chi2_full, bounds, args=(ener_f, ic, r_sun),
#                                       strategy='best2bin', maxiter=300, popsize=300, tol=5.0e-8,
#                                       atol=0.0, disp=True, polish=True, workers=-1)
# param_fitted = opt.x
# np.savetxt(param_file, param_fitted, delimiter=',')
# w_0 = param_fitted
w_0 = np.loadtxt(param_file)


theta_0 = w_0[0]
d_theta = w_0[1]
beta_0 = w_0[2]
W_0 = theta_0 + d_theta
pot_list = pot_model(ener_f, theta_0, W_0, beta_0)
phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel, x, y, z, v_circ = orbit_model(ic[0], ic[1], ic[2], ic[3], ic[4],
                                                                         ic[5], pot_list, r_sun)

print('Model parameters:')
print('theta_0, W_0, beta_0, d_theta =', theta_0, W_0, beta_0, d_theta)
print('IC:', ic)
chi2_stream(w_0, ener_f, ic, r_sun)
chi2_full(w_0, ener_f, ic, r_sun)

# Plot in galactocentric coordinates
fig = plt.figure(figsize=(10, 10))

plt.scatter(x, y, s=0.1, marker='o', color='red')
plt.scatter(x, z, s=0.1, marker='o', color='blue')
plt.xlim(-15, 10)
plt.ylim(-15, 20)
plt.grid()
plt.legend()
plt.tight_layout()
fig.savefig("plots/orbit_fit_pot_from_IbataPolysGaiaDR2_chi2full.png")
plt.show()

# Plots in the sky using the GD-1 frame
fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, sharex=True, figsize=(10, 50))


# Sky position
ax1.set_title('Ibata+20 (Polynomials) data fitted by RAR halo plus barions')
ax1.scatter(phi_1.wrap_at(180*u.deg), phi_2, s=0.1, marker='o', color='red', label='RAR+barions fit')
ax1.plot(Iba_sky['phi_1'], Iba_sky['phi_2'], color='blue', label='GD-1 (Ibata+2020)')
ax1.set_ylim(-4, 2)
ax1.set_ylabel(r'$\phi_2$ [degrees]')
ax1.legend()

# Heliocentric distance
ax2.scatter(phi_1, d_hel, s=0.1, marker='o', color='red')
ax2.plot(Iba_sky['phi_1'], Iba_sky['d_hel'], color='blue', label='Stream\n(Ibata+2020)')
ax2.plot(Iba_sky['phi_1'], r_hel, color='green', label='Integrated from v_hel')
ax2.set_ylim(7, 12)
ax2.set_ylabel(r'$D$ [kpc]')

# Proper motion along RA
ax3.scatter(phi_1, mu_ra, s=0.5, color='red')
ax3.plot(Iba_sky['phi_1'], Iba_sky['mu_ra'], color='blue', label='Stream\n(Ibata+2020)')
ax3.set_ylabel(r'$\mu_\alpha$ [mas yr$^{-1}$]')
ax3.set_ylim(-9, 0)


# Proper motion along DEC
ax4.scatter(phi_1, mu_dec, s=0.5, color='red', label='RAR+barions fit')
ax4.plot(Iba_sky['phi_1'], Iba_sky['mu_dec'], color='blue')
ax4.set_ylabel(r'$\mu_\delta$ [mas yr$^{-1}$]')

# Heliocentric radial velocity
ax5.scatter(phi_1.wrap_at(180*u.deg), v_hel, s=0.1, marker='o', color='red')
ax5.plot(Iba_sky['phi_1'], Iba_sky['v_hel'], color='blue', label='Stream\n(Ibata+2020)')
ax5.set_ylim(-300, 300)
ax5.set_ylabel(r'$v_{\rm{LOS}}$ [km s$^{-1}$]')

# Galactocentric distance
d_gal = np.sqrt(x*x+y*y+z*z)
ax6.scatter(phi_1.wrap_at(180*u.deg), d_gal, s=0.1, marker='o', color='red')
ax6.plot(Iba_sky['phi_1'], Iba_d_gal, color='blue', label='Stream\n(Ibata+2020)')
ax6.plot(Iba_sky['phi_1'], Iba_r_gal, color='green', label='Integrated from v_hel')
ax6.set_ylim(13, 17)
ax6.set_ylabel(r'$d_{\rm{Gal}}$ [kpc]')


plt.xlabel(r'$\phi_1$ [degrees]')
plt.xlim(IbaPoly.limit[0], IbaPoly.limit[1])
fig.savefig("plots/sky_fit_pot_from_IbataPolysGaiaDR2_chi2full.png")
plt.show()
