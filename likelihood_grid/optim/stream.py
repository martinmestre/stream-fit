"""
Likelihood.

Computation of the orbit model for RAR + barionic component.
Using Ibata+20 polinomial data-fits.

Author: Martín Mestre.
"""


import numpy as np
import astropy.coordinates as coord
import pandas as pd
from astropy import units as u
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy import optimize
from scipy.signal import argrelextrema
import GD1Koposov10_class as GD1_class
import potential_classes as pot
import config as cfg


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


def v_circ_Newton(r, halo):
    """Circular velocity for the rar model."""
    G_u = 4.3009e-6  # kpc (km/s)^2 M_sun^-1
    return np.sqrt(G_u*halo.mass_wrap(r)/r)


def get_core(halo):
    """Core mass."""
    arg_max = argrelextrema(v_circ_Newton(halo.r_s, halo), np.greater)
    arg_first_max = arg_max[0][0]
    r_cand = halo.r_s[arg_first_max]
    bounds = np.array([r_cand*0.5, r_cand*1.5])
    r_core = optimize.fminbound(lambda r: -v_circ_Newton(r, halo),
                                bounds[0], bounds[1], xtol=0.5e-12, maxfun=1000)
    m_core = halo.mass_wrap(r_core)
    return r_core, m_core


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
    """Orbit definition."""
    """
    The orbit_model here defined works with sky coordinates
    at input and sky-cartesian at output.
    """

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
    gd1_coord = galac_coord.transform_to(GD1_class.GD1Koposov10())
    phi_1 = gd1_coord.phi1
    phi_2 = gd1_coord.phi2
    d_hel = gd1_coord.distance
    v_hel = gd1_coord.radial_velocity
    # Unused block of code
    # mu_phi_1 = gd1_coord.pm_phi1_cosphi2/np.cos(phi_2)  #not used by Ibata
    # mu_phi_2 = gd1_coord.pm_phi2
    # return phi_1, phi_2, d_hel, v_hel, mu_phi_1, mu_phi_2
    # Transformation to ICRS coordinates
    icrs_coord = galac_coord.transform_to(coord.ICRS())
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

# Core constraint
m_core_const = 3.5e6  # M_sun


def chi2_stream(theta_0, d_theta, beta_0, ener_f, ic, r_sun):
    """Chi^2 stream function."""
    import wrap
    W_0 = theta_0 + d_theta
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

    print(theta_0, d_theta, beta_0, ener_f, '--- chi2_stream = ', np.sum(sum))
    print("::: r_core = ", r_core, " ::: m_core = ", mass_core/1.e6, "x10⁶ M_sun")
    return np.sum(sum)


def chi2_core(theta_0, d_theta, beta_0, ener_f):
    """Chi^2 core function ."""
    W_0 = theta_0 + d_theta
    param = np.array([ener_f, theta_0, W_0, beta_0])
    halo = pot.RAR(param)
    r_core, mass_core = get_core_GR(halo)
    chisq_core = (mass_core - m_core_const)**2/(0.01*m_core_const)**2
    print(theta_0, d_theta, beta_0, ener_f, '--- chi2_core = ', chisq_core)
    print("::: r_core = ", r_core, " ::: m_core = ", mass_core/1.e6, "x10⁶ M_sun")
    return chisq_core


def chi2_full(theta_0, d_theta, beta_0, ener_f, ic, r_sun):
    """Chi^2 full (stream+core) function."""
    import wrap
    W_0 = theta_0 + d_theta
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
    return np.sum(sum)


def chi2_stream_potlist(pot_list, ic, r_sun):
    """Chi^2 stream function using the potential object as argument."""
    import wrap
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

    print("chi2_stream_pot_list = ", np.sum(sum), "  r_core = ", r_core,
          "  m_core = ", mass_core/1.e6, "x10⁶ M_sun")
    return np.sum(sum)
