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


def pot_model(ener_f, theta_0, W_0, beta_0, alpha, scale, scaleb):
    """Potential model including barions and dark matter halo."""
    # Set barionic potentials
    M_gal = 2.32e7  # Msun
    bulge = pot.Plummer(460.0*M_gal, 0.3)
    thin = pot.MiyamotoNagai(1700.0*M_gal*alpha, 5.3*scale, 0.25*scaleb)
    thick = pot.MiyamotoNagai(1700.0*M_gal*alpha, 2.6*scale, 0.8*scaleb)
    # Set halo potential
    param = np.array([ener_f, theta_0, W_0, beta_0])
    halo = pot.RAR(param)
    return bulge, thin, thick, halo


def orbit_model(alpha, delta, distance, mu_alpha, mu_delta, v_los, pot_list):
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
    mu_ra = icrs_coord.pm_ra_cosdec  # / np.cos(icrs_coord.dec) # It seems that Ibata gives pm_ra_cosdec
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


def chi2(w_0, ener_f, beta_0, ic):
    """Chi^2 function."""
    import wrap
    theta_0 = w_0[0]
    d_theta = w_0[1]
    alpha = w_0[2]
    scale = w_0[3]
    scaleb = w_0[4]
    W_0 = theta_0 + d_theta
    pot_list = pot_model(ener_f, theta_0, W_0, beta_0, alpha, scale, scaleb)
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

    print('chi^2 =', np.sum(sum), ' v_circ = ', v_circ, ' theta_0=', theta_0, ' W_0=', W_0,
          ' alpha=', alpha, ' scale=', scale, 'scaleb=', scaleb)
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

# Taking the initial conditions from the Galpy fit with fixed MW2014 potential.
ic = np.array([1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
              -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01])


# Parameters
param_file = 'param_fit_pot-slice-barions_from_IbataPolysGaiaDR2-data.txt'
r_sun = 8.0    # 8.122
ener_f = 56.0  # keV
# d_theta = 28.5751
beta_0 = 1.25e-5  # 1.1977e-5


# Optimization
# bounds = ((34, 39), (25, 30), (0.85, 1.15), (0.9, 1.1), (0.9, 1.1))
# opt = optimize.differential_evolution(chi2, bounds, args=(ener_f, beta_0, ic),
#                                      strategy='best2bin', maxiter=40, popsize=40, tol=5.0e-6,
#                                      atol=0.5e-6, disp=True, polish=True, workers=-1)
# param_fitted = opt.x
# np.savetxt(param_file, param_fitted, delimiter=',')
# w_0 = param_fitted
w_0 = np.loadtxt(param_file)
chi2(w_0, ener_f, beta_0, ic)

theta_0 = w_0[0]
d_theta = w_0[1]
alpha = w_0[2]
scale = w_0[3]
scaleb = w_0[4]
W_0 = theta_0 + d_theta
pot_list = pot_model(ener_f, theta_0, W_0, beta_0, alpha, scale, scaleb)
phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel, x, y, z, v_circ = orbit_model(ic[0], ic[1], ic[2], ic[3], ic[4],
                                                                         ic[5], pot_list)

print('Model parameters:')
print('theta_0, W_0, beta_0, alpha =', theta_0, W_0, beta_0, alpha)
print('IC:', ic)

# Plot in galactocentric coordinates
fig = plt.figure(figsize=(10, 10))

plt.scatter(x, y, s=0.1, marker='o', color='red')
plt.scatter(x, z, s=0.1, marker='o', color='blue')
plt.xlim(-15, 10)
plt.ylim(-15, 20)
plt.grid()
plt.legend()
plt.tight_layout()
# plt.show()
fig.savefig("plots/orbit_pot-slice-barions_from_IbataPolysGaiaDR2-data.png")


# Plots in the sky using the GD-1 frame
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex=True, figsize=(7, 35))


# Sky position
ax1.set_title('Ibata+20 data with RAR and barions')
ax1.scatter(phi_1.wrap_at(180*u.deg), phi_2, s=0.1, marker='o', color='red')
ax1.plot(Iba_sky['phi_1'], Iba_sky['phi_2'], color='blue', label='Stream\n(Ibata+2020)')
ax1.set_ylim(-4, 2)
ax1.set_ylabel(r'$\phi_2$ [degrees]')


# Heliocentric radial velocity
ax2.scatter(phi_1.wrap_at(180*u.deg), v_hel, s=0.1, marker='o', color='red')
ax2.plot(Iba_sky['phi_1'], Iba_sky['v_hel'], color='blue', label='Stream\n(Ibata+2020)')
ax2.set_ylim(-300, 300)
ax2.set_ylabel(r'$v_{\rm{LOS}}$ [km s$^{-1}$]')

# Heliocentric distance
ax3.scatter(phi_1, d_hel, s=0.1, marker='o', color='red')
ax3.plot(Iba_sky['phi_1'], Iba_sky['d_hel'], color='blue', label='Stream\n(Ibata+2020)')
ax3.set_ylim(7, 12)
ax3.set_ylabel(r'$D$ [kpc]')

# Proper motion along RA
ax4.scatter(phi_1, mu_ra, s=0.5, color='red')
ax4.plot(Iba_sky['phi_1'], Iba_sky['mu_ra'], color='blue', label='Stream\n(Ibata+2020)')
ax4.set_ylabel(r'$\mu_\alpha$ [mas yr$^{-1}$]')
ax4.set_ylim(-20, 0)
ax4.legend()

# Proper motion along DEC
ax5.scatter(phi_1, mu_dec, s=0.5, color='red')
ax5.plot(Iba_sky['phi_1'], Iba_sky['mu_dec'], color='blue', label='Stream\n(Ibata+2020)')
ax5.set_ylabel(r'$\mu_\delta$ [mas yr$^{-1}$]')
ax5.legend()

plt.xlabel(r'$\phi_1$ [degrees]')
plt.xlim(IbaPoly.limit[0], IbaPoly.limit[1])
# plt.tight_layout()

plt.show()
fig.savefig("plots/sky_pot-slice-barions_from_IbataPolysGaiaDR2-data.png")
