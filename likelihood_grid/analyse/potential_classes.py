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
        cfg.r_max = np.amax(self.r)
        cfg.mass_spline = InterpolatedUnivariateSpline(self.r, self.mass, k=4)
        cfg.nu_spline = InterpolatedUnivariateSpline(self.r, self.nu, k=4)

    def mass_wrap(self, r):
        """Wrap."""
        if(r < cfg.r_max):
            return cfg.mass_spline(r)
        else:
            return cfg.mass_spline(cfg.r_max)

    def accel(self, x, y, z):
        """Acceleration."""
        r = np.sqrt(x*x+y*y+z*z)
        fac = -g*self.mass_wrap(r)/r**3
        return fac*np.array([x, y, z])

    def dnu(self, r):
        """Derive nu."""
        deriv = cfg.nu_spline.derivative(1)
        return deriv(r)
