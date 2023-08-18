"""RAR model.

Metric convention:
g_00 = e^(nu)
g_11 = -e^(lambda)

"""

from scipy.integrate import solve_ivp
import numpy as np
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline


def model(param):
    """RAR model as a function."""
    def fermi(eps, alpha_r, beta_r, eps_r):
        """
        Fermi distribution.

        It also depends on alpha(r), beta(r) and eps(r).
        No explicit dependence on r.
        """
        up = 1.0 - np.exp((eps-eps_r)/beta_r)
        down = 1.0 + np.exp((eps-alpha_r)/beta_r)
        return up/down

    def g_rho(eps, alpha_r, beta_r, eps_r):
        """Integrand for the density."""
        return eps*eps*np.sqrt((eps-1.0)*(eps+1.0))*fermi(eps, alpha_r, beta_r, eps_r)

    def g_P(eps, alpha_r, beta_r, eps_r):
        """Integrand for the pressure."""
        return np.sqrt((eps-1.0)*(eps+1.0))**3*fermi(eps, alpha_r, beta_r, eps_r)

    def density(n_step, alpha_r, beta_r, eps_r):
        """Density."""
        if(eps_r > 1.0 + float(n_step)*machine_eps):
            eps = np.linspace(1.0, eps_r, n_step)
            return a * integrate.simps(g_rho(eps, alpha_r, beta_r, eps_r), eps, even='avg')
        else:
            return 0.0

    def pressure(n_step, alpha_r, beta_r, eps_r):
        """Pressure."""
        if(eps_r > 1.0 + float(n_step)*machine_eps):
            eps = np.linspace(1.0, eps_r, n_step)
            return b * integrate.simps(g_P(eps, alpha_r, beta_r, eps_r), eps, even='avg')
        else:
            return 0.0

    def border_energy_cutoff(t, u):
        """Border event function: energy_cutoff. Not used anymore."""
        exponential = np.exp(-0.5*(u[1]-nu_0))
        eps_r = eps_0*exponential
        return eps_r - 1.0

    border_energy_cutoff.terminal = True
    border_energy_cutoff.direction = -1

    def border_density(t, u):
        """Border event function: density."""
        exponential = np.exp(-0.5*(u[1]-nu_0))
        alpha_r = alpha_0*exponential
        beta_r = beta_0*exponential
        eps_r = eps_0*exponential
        rho = density(n_eos, alpha_r, beta_r, eps_r)
        return rho*to_astro - 1.0e-10  # M_sun/pc^3

    border_density.terminal = True
    border_density.direction = -1

    def dnu_dt(t, u):
        """Derive nu respect to t."""
        # Computing density and pressure:
        exponential = np.exp(-0.5*(u[1]-nu_origin))
        alpha_r = alpha_0*exponential
        beta_r = beta_0*exponential
        eps_r = eps_0*exponential
        if(eps_r < 1.0):           # This is important to avoid NaNs and to define the halo border.
            eps_r = 1.0

        # rho = density(n_eos, alpha_r, beta_r, eps_r)
        P = pressure(n_eos, alpha_r, beta_r, eps_r)

        d_nu = (np.exp(u[0]) + np.exp(2.0*t)*P/(rho_rel*c*c))/(1.0-np.exp(u[0]))
        return d_nu

    def density_astro_kpc(nu):
        """Density to ."""
        exponential = np.exp(-0.5*(nu-nu_0))
        alpha_r = alpha_0*exponential
        beta_r = beta_0*exponential
        eps_r = eps_0*exponential
        return density(n_eos, alpha_r, beta_r, eps_r)*to_astro_kpc

    def tov(t, u):
        """
        Tolman-Oppenheimer-Volkoff equation.

        t: independent variable, t=ln(r/R).
        u: 5D array whose components are:
        u[0] = z          with z=ln(Psi) and Psi = (M(r)/M)*(R/r).
        u[1] = nu            = metric potential.
        where:
        T = temperature(r).
        theta = degeneracy parameter = mu /(k*T).
        mu = chemical potential with rest mass subtracted off.
        W = cutoff parameter = E_c / (k*T).
        E_c = cutoff energy.
        M = mass scaling factor.
        R = radius scaling factor.
        """
        # Computing density and pressure:
        exponential = np.exp(-0.5*(u[1]-nu_0))
        alpha_r = alpha_0*exponential
        beta_r = beta_0*exponential
        eps_r = eps_0*exponential
        if(eps_r < 1.0):           # This is important to avoid NaNs and to define the halo border.
            eps_r = 1.0

        rho = density(n_eos, alpha_r, beta_r, eps_r)
        P = pressure(n_eos, alpha_r, beta_r, eps_r)

        d_z = np.exp(2.0*t-u[0])*rho/rho_rel - 1.0
        d_nu = (np.exp(u[0]) + np.exp(2.0*t)*P/(rho_rel*c*c))/(1.0-np.exp(u[0]))

        return [d_z, d_nu]

    # Set constants
    kpc2pc = 1000.0
    cm2kpc = 1.0/3.08567758e21
    cm2pc = cm2kpc*kpc2pc
    g2Msun = 1.0/1.98847e33
    to_astro = g2Msun/(cm2pc**3)
    to_astro_kpc = g2Msun/(cm2kpc**3)
    pi = np.pi
    c = 2.99792458e+10
    G = 6.67430e-8
    h = 6.6260755e-27
    g = 2.0
    m_e = 9.1093837015e-28  # cgs
    ener_e = 510.99895  # keV
    ener_f = param[0]  # 48.0  keV
    m = (ener_f/ener_e)*m_e  # cgs
    rho_rel = (g*m**4/h**3)*(pi*c*c)**1.5
    R = c/np.sqrt(8.0*pi*G*rho_rel)
    M = 4.0*pi*R**3*rho_rel
    a = 4.0*rho_rel/np.sqrt(pi)
    b = a*c*c/3.0
    n_eos = 2**10+1
    tau = 1.0e-15
    min_r = 1.0e-16  # kpc
    max_r = 1.0e3  # kpc
    machine_eps = np.finfo(float).eps

    # Set initial conditions
    psi_0 = 2.0*tau
    z_0 = np.log(psi_0)
    nu_0 = 2.0*tau
    theta_0 = param[1]  # 35.55
    W_0 = param[2]  # 80.0
    beta_0 = param[3]  # 1.0e-5
    alpha_0 = 1.0+beta_0*theta_0
    eps_0 = 1.0 + beta_0*W_0
    rho_0 = density(n_eos, alpha_0, beta_0, eps_0)
    u_0 = [z_0, nu_0]
    t_0 = np.log(np.sqrt(6.0*tau*rho_rel/rho_0))
    t_f = np.log(max_r/cm2kpc/R)

    # Solving the TOV system
    sol = solve_ivp(tov, [t_0, t_f], u_0, method='LSODA',
                    events=(border_density),
                    rtol=5.0e-14, atol=0.0)

    # Defining physical variables
    r = np.exp(sol.t)*R
    z = sol.y[0]
    nu = sol.y[1]
    psi = np.exp(z)
    mass = psi*M/R*r
    exponential = np.exp(-0.5*(nu-nu_0))
    beta = beta_0*exponential
    rho = np.array([density_astro_kpc(nu[i]) for i in range(0, len(sol.t))])


    # Shift the metric:
    nu_origin = 2.0*np.log(np.sqrt(1.0-psi[-1])*beta[-1]/beta[0])
    nu = nu + nu_origin*np.ones(len(nu))


    # In astrophysical units
    r = r*cm2kpc  # kpc
    mass = mass*g2Msun  # M_sun
    bool_r = (r > min_r)
    r = r[bool_r]
    mass = mass[bool_r]
    nu = nu[bool_r]
    r_t = sol.t[bool_r]
    z = z[bool_r]
    rho = rho[bool_r]

    dnu_dr = np.array([dnu_dt(r_t[i], [z[i], nu[i]])/r[i] for i in range(0, len(r))])

    # Uncomment only if we need to check the Density.
    # mass_spline = InterpolatedUnivariateSpline(r, mass, k=4)  # Allows easy computation of derivatives
    # def rho_spline(r):
    #     """Density function, afterwards."""
    #     deriv = mass_spline.derivative(1)
    #     return deriv(r)/(4.0*np.pi*r*r)
    # rho = rho_spline(r)

    return r, mass, nu, dnu_dr, rho
