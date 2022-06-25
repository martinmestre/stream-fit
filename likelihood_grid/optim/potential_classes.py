"""Potential classes."""

# Gravitational Potential Classes.
import config as self
import numpy as np
import model_def
from scipy.interpolate import InterpolatedUnivariateSpline

g = self.g


class MiyamotoNagai:
    """Miyamoto-Nagai disk."""

    def __init__(self, M, a, b):
        """Init."""
        self.M = M
        self.a = a
        self.b = b
        return None

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
        r, mass, nu = model_def.model(self.param)
        isnan = np.argwhere(np.isnan(mass))
        if (np.any(isnan)):
            k = np.argwhere(np.isnan(mass))[0][0]
        else:
            k = -1
        self.r = r[0:k]
        self.mass = mass[0:k]
        self.nu = nu[0:k]
        self.r_max = np.amax(self.r)
        self.mass_spline = InterpolatedUnivariateSpline(self.r, self.mass, k=4)
        self.nu_spline = InterpolatedUnivariateSpline(self.r, self.nu, k=4)
        self.dnu_spline = self.nu_spline.derivative(1)

    def mass_wrap(self, r):
        """Wrap."""
        if(r < self.r_max):
            return self.mass_spline(r)
        else:
            return self.mass_spline(self.r_max)

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
        if(r < self.r_max):
            return self.dnu_spline(r)
        else:
            G_u = 4.3009e-6  # kpc (km/s)^2 M_sun^-1
            c = 2.997925e5  # km s^-1
            m_total = self.mass_spline(self.r_max)
            r_Schwar = 2.0*G_u*m_total/c**2
            dnu_Schwar = r_Schwar*np.exp(-self.nu_wrap(r))/r**2
            return dnu_Schwar

    def accel(self, x, y, z):
        """Acceleration."""
        r = np.sqrt(x*x+y*y+z*z)
        fac = -g*self.mass_wrap(r)/r**3
        return fac*np.array([x, y, z])
