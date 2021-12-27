"""RAR model."""

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.optimize import least_squares
import astropy.constants as cons
import astropy.units as u
from fractions import Fraction


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

    # Input arguments
    ener, theta_0, W_0, beta_0 = param

    # Constants
    pi = np.pi
    g = 2.0
    to_astro = (u.g).to(u.Msun)/((u.cm).to(u.pc))**3
    c = cons.c.cgs   # 2.99792458e+10
    G = cons.G.cgs   # 6.67430e-8
    h = cons.h.cgs   # 6.6260755e-27
    m = (ener/c**2).cgs
    rho_rel = ((g*m**4/h**3)*(pi*c*c)**Fraction(3, 2)).cgs
    R = c/np.sqrt(8.0*pi*G*rho_rel)
    M = 4.0*pi*R**3*rho_rel
    a = (4.0*rho_rel/np.sqrt(pi)).value
    b = (a*c*c).value
    n_eos = 2**10+1
    tau = 1.0e-15
    min_r = 1.0e-16*u.kpc
    max_r = 1.0e3*u.kpc
    machine_eps = np.finfo(float).eps

    # Set initial conditions
    psi_0 = 2.0*tau
    z_0 = np.log(psi_0)
    nu_0 = tau
    alpha_0 = 1.0+beta_0*theta_0
    eps_0 = 1.0 + beta_0*W_0
    rho_0 = density(n_eos, alpha_0, beta_0, eps_0)*u.g/u.cm**3

    u_0 = [z_0, nu_0]
    t_0 = np.log(np.sqrt(6.0*tau*rho_0/rho_rel))
    t_f = np.log(max_r/R)

    c = c.value
    rho_rel = rho_rel.value

    # Solving the TOV system in adimensional units
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

    # Shift the metric:
    nu_origin = 2.0*np.log(np.sqrt(1.0-psi[-1])*beta[-1]/beta[0])
    nu = nu + nu_origin*np.ones(len(nu))

    # In astrophysical units
    r = r.to(u.kpc)  # kpc
    mass = mass.to(u.Msun)  # M_sun
    bool_r = (r > min_r)
    r = r[bool_r]
    mass = mass[bool_r]
    nu = nu[bool_r]

    return r, mass, nu
