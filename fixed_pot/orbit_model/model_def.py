# Fit RAR model to mass(r) data points.
#--------------------------------------


from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import least_squares


# Model definition
#-----------------


def fermi(eps,alpha_r,beta_r,eps_r):                # it also depends on alpha(r), beta(r) and eps(r).
                           # no explicit dependence on r.
    up = 1.0 - np.exp((eps-eps_r)/beta_r)
    down = 1.0 + np.exp((eps-alpha_r)/beta_r)
    return up/down

def g_rho(eps,alpha_r,beta_r,eps_r):
    return eps*eps*np.sqrt(eps*eps-1.0)*fermi(eps,alpha_r,beta_r,eps_r)

def g_P(eps,alpha_r,beta_r,eps_r):
    return np.sqrt(eps*eps-1.0)**3*fermi(eps,alpha_r,beta_r,eps_r)

def density(n_step,alpha_r,beta_r,eps_r):
    eps=np.linspace(1.0,eps_r,n_step)
    return a * integrate.simps(g_rho(eps,alpha_r,beta_r,eps_r), eps, even='avg')

def pressure(n_step,alpha_r,beta_r,eps_r):
    eps=np.linspace(1.0,eps_r,n_step)
    return b * integrate.simps(g_P(eps,alpha_r,beta_r,eps_r), eps, even='avg')

def tov(t,u):
    # t: independent variable, t=ln(r/R).
    # u: 5D array whose components are:
    # u[0] = z          with z=ln(\Psi) and \Psi = (M(r)/M)*(R/r)
    # u[1] = \nu            = metric potential
    # where: 
    # T = temperature(r)
    # \theta = degeneracy param. = \mu /(k*T)
    # \mu = chemical potential with rest mass subtracted off 
    # W = cutoff param. = E_c / (k*T)
    # E_c = cutoff energy
    # M = mass scaling factor
    # R = radius scaling factor
    
    
    # Computing density and pressure: 
    exponential =np.exp(-0.5*(u[1]-nu_0))
    alpha_r = alpha_0*exponential
    beta_r = beta_0*exponential
    eps_r= eps_0*exponential
    
    rho = density(n_eos,alpha_r,beta_r,eps_r)
    P   = pressure(n_eos,alpha_r,beta_r,eps_r)
    
    d_z = np.exp(2.0*t-u[0])*rho/rho_rel - 1.0
    d_nu  = ( np.exp(u[0]) + np.exp(2.0*t)*P/(rho_rel*c*c) )/(1.0-np.exp(u[0]))
        
    
    return [d_z,d_nu]


def model(param):
    global nu_0,alpha_0,beta_0,eps_0
    global a,b,c
    global n_eos,rho_rel
    global R,M
    
    # Set constants
    cm2kpc=1.0/3.08567758e21
    g2Msun=1.0/1.98847e33
    pi = np.pi
    c = 2.99792458e+10
    G = 6.67430e-8
    h = 6.6260755e-27
    g = 2.0
    m_e = 9.1093837015e-28  # cgs
    ener_e = 510.99895 # keV
    ener_f = param[0] #     48.0  keV
    m = (ener_f/ener_e)*m_e   # cgs
    rho_rel = (g*m**4/h**3)*(pi*c*c)**1.5
    R = c/ np.sqrt(8.0*pi*G*rho_rel)
    M = 4.0*pi*R**3*rho_rel
    a = 4.0*rho_rel/np.sqrt(pi)
    b = a*c*c
    n_eos = 10**4
    n_tov =  10**6
    tau = 1.0e-16


    # Set initial conditions
    psi_0 = 2.0*tau
    z_0=np.log(psi_0)
    nu_0 = tau
    theta_0= param[1]  # 35.55
    W_0 = param[2] #80.0
    beta_0 = param[3]# 1.0e-5
    alpha_0 =1.0+beta_0*theta_0
    eps_0 = 1.0+ beta_0*W_0
    rho_0 = density(n_eos,alpha_0,beta_0,eps_0)

    u_0 = [z_0,nu_0]
    t_0 = np.log(np.sqrt(6.0*tau*rho_0/rho_rel))   
    t_f = np.log(1.0e+11)
    t_eval = np.linspace(t_0,t_f,n_tov) 

    print('r_0=', np.exp(t_0)*R*cm2kpc)
    print('r_f=', np.exp(t_f)*R*cm2kpc)
    print('t_0=',t_0)
    print('R=',R)

    print(theta_0,W_0,beta_0)
  
    # Solving the TOV system
    sol = solve_ivp(tov, [t_0, t_f], u_0, t_eval=t_eval,method='LSODA',rtol=5.0e-14,atol=0.5e-14)

    # Defining physical variables
    r= np.exp(sol.t)*R
    z=sol.y[0]
    nu = sol.y[1]

    psi=np.exp(z)
    mass= psi*M/R*r
    exponential = np.exp(nu-nu_0)
    alpha = alpha_0*exponential
    beta = beta_0*exponential
    eps = eps_0*exponential

#    plt.plot(r/R,mass/M)
#    plt.plot(r[:-1]/R,np.diff(mass)/M/np.diff(r)/(4.0*pi*r[:-1]**2))
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.show()
 
    return r*cm2kpc,mass*g2Msun







