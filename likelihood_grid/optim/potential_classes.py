"""Potential classes."""

# Gravitational Potential Classes.
import config as cfg
import numpy as np
import model_def
from scipy.interpolate import InterpolatedUnivariateSpline

g = cfg.g


class MiyamotoNagai:
    """Miyamoto-Nagai disk."""

    def __init__(self, M, a, b):
        """Init."""
        self.M = M
        self.a = a
        self.b = b

    def accel(self, x, y, z):
        """Acceleration."""
        fac = -g*self.M / (x*x + y*y + (self.a + np.sqrt(self.b**2 + z*z))**2)**(1.5)
        return np.array([fac*x, fac*y, fac*z*(1.0+self.a/np.sqrt(self.b**2+z*z))])


class Plummer:
    """Plummer sphere."""

    def __init__(self, M, b):
        """Init."""
        self.M = M
        self.b = b

    def accel(self, x, y, z):
        """Acceleration."""
        fac = -g*self.M/(self.b**2+x*x+y*y+z*z)**(1.5)
        return np.array([fac*x, fac*y, fac*z])


class RAR:
    """RAR halo."""

    def __init__(self, param):
        """Init."""
        self.param = param
        r, mass, nu, dnu = model_def.model(self.param)
        isnan = np.argwhere(np.isnan(mass))
        if (np.any(isnan)):
            k = np.argwhere(np.isnan(mass))[0][0]
        else:
            k = -1
        r_s = r[0:k]
        mass_s = mass[0:k]
        nu_s = nu[0:k]
        print('radio=', r_s[-1], ' kpc')
        print('masa=', mass_s[-1]/1.e11, ' x10^11 solar masses')
        self.r_max = np.amax(r_s)
        self.r_s = r_s
        self.mass_spline = InterpolatedUnivariateSpline(r_s, mass_s, k=4)
        self.nu_spline = InterpolatedUnivariateSpline(r_s, nu_s, k=4)
        # self.dnu_spline = self.nu_spline.derivative(1)
        self.dnu_integral_spl = InterpolatedUnivariateSpline(r_s, dnu, k=4)

    def mass_wrap(self, r):
        """Wrap."""
        if(np.isscalar(r)):
            if(r < self.r_max):
                mass = self.mass_spline(r)
            else:
                mass = self.mass_spline(self.r_max)
        else:
            mass = np.zeros(len(r))
            for i in range(0, len(r)):
                if(r[i] < self.r_max):
                    mass[i] = self.mass_spline(r[i])
                else:
                    mass[i] = self.mass_spline(self.r_max)
        return mass

    def nu_wrap(self, r):
        """Metric nu wrapped by Schwarschild."""
        if(r < self.r_max):
            return self.nu_spline(r)
        else:
            G_u = 4.3009e-6  # kpc (km/s)^2 M_sun^-1
            c = 2.997925e5  # km s^-1
            m_total = self.mass_spline(self.r_max)
            r_Schwar = 2.0*G_u*m_total/c**2
            return np.log(1.0-r_Schwar/r)

    def dnu_wrap(self, r):
        """Wrap dnu with Schwarschild metric."""
        if(np.isscalar(r)):
            if(r < self.r_max):
                dnu = self.dnu_integral_spl(r)
            else:
                G_u = 4.3009e-6  # kpc (km/s)^2 M_sun^-1
                c = 2.997925e5  # km s^-1
                m_total = self.mass_spline(self.r_max)
                r_Schwar = 2.0*G_u*m_total/c**2
                dnu_Schwar = r_Schwar*np.exp(-self.nu_wrap(r))/r**2
                dnu = dnu_Schwar
        else:
            dnu = np.zeros(len(r))
            for i in range(0, len(r)):
                if(r[i] < self.r_max):
                    dnu[i] = self.dnu_integral_spl(r[i])
                else:
                    G_u = 4.3009e-6  # kpc (km/s)^2 M_sun^-1
                    c = 2.997925e5  # km s^-1
                    m_total = self.mass_spline(self.r_max)
                    r_Schwar = 2.0*G_u*m_total/c**2
                    dnu_Schwar = r_Schwar*np.exp(-self.nu_wrap(r[i]))/r[i]**2
                    dnu[i] = dnu_Schwar
        return dnu

    def accel(self, x, y, z):
        """Acceleration."""
        r = np.sqrt(x*x+y*y+z*z)
        fac = -g*self.mass_wrap(r)/r**3
        return fac*np.array([x, y, z])
