# Gravitational Potential Classes.
import config as cfg
import numpy as np
import model_def 
from scipy.interpolate import InterpolatedUnivariateSpline

g = cfg.g


class MiyamotoNagai:
	def __init__(self,M,a,b):
		self.M = M
		self.a = a
		self.b = b
	def accel(self,x,y,z):
		fac = -g*self.M / (x*x + y*y + (self.a + np.sqrt(self.b**2 + z*z))**2)**(1.5)
		return np.array([fac*x, fac*y, fac*z*(1.0+self.a/np.sqrt(self.b**2+z*z)) ])


class Plummer:
	def __init__(self,M,b):
		self.M = M
		self.b = b
	def accel(self,x,y,z):
		fac = -g*self.M/(self.b**2 +x*x+y*y+z*z )**(1.5)
		return np.array([fac*x, fac*y, fac*z])


class RAR:
	def __init__(self,param):
		self.param = param
		r,mass,nu = model_def.model(self.param)
		isnan = np.argwhere(np.isnan(mass))
		if (np.any(isnan)):
			k = np.argwhere(np.isnan(mass))[0][0]
		else:
			k=-1
		r_s=r[0:k]
		mass_s=mass[0:k]
		cfg.r_max = np.amax(r_s)
		cfg.mass_spline = InterpolatedUnivariateSpline(r_s,mass_s,k=4)

	def mass_wrap(self,r):
		if(r<cfg.r_max):
			return cfg.mass_spline(r)
		else:
			return cfg.mass_spline(cfg.r_max)

	def accel(self,x,y,z):
		r=np.sqrt(x*x+y*y+z*z)
		fac= -g*self.mass_wrap(r)/r**3
		return fac*np.array([x,y,z])
    
    
