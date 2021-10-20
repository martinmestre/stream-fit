
# Gravitational accelerations
#-----------------------------
import scipy.special as spec
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import config as cfg

# Parameters
g= cfg.g


# Milky Way
#----------



# RAR (halo)
def rho_rar_mw(r):
    if(r<cfg.r_max_mw):
        return rho_spline(r)
    else:
        return rho_spline(cfg.r_max_mw)

def mass_rar_mw(r):
    if(r<cfg.r_max_mw):
        return cfg.mass_spline(r)
    else:
        return cfg.mass_spline(cfg.r_max_mw)
    
def accel_rar_mw(x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    fac= -g*mass_rar_mw(r)/r**3
    return fac*np.array([x,y,z])


# NFW (halo)
def f(x):
    return np.log(1.0+x)- x/(1.0+x)

def rho_nfw_mw(r):
    rho_0 = cfg.m_200/(4.0*np.pi*f(c)*cfg.a_nfw**3)
    return rho_0/(r/cfg.a_nfw)/(1+r/cfg.a_nfw)**2

def mass_nfw_mw(r):
    return (cfg.m_200/f(c))*f(r/cfg.a_nfw)


def accel_nfw_mw(x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    fac= -g*mass_nfw_mw(r)/r**3
    return fac*np.array([x,y,z])


# Logarithmic (halo)
def accel_logarithmic(x,y,z):
    r2 = x*x+y*y+z*z/(cfg.q_phi**2)
    return -cfg.v_c**2/r2*np.array([x,y,z/(cfg.q_phi**2)])


# Plummer sphere (bulge)
def accel_plum_bulge(x,y,z):
    fac = -g*cfg.M_plu/(cfg.r_plu**2 +x*x+y*y+z*z )**(1.5)
    return np.array([fac*x, fac*y, fac*z])

# Miyamoto-Nagai (disk)
def accel_mn_disk(x,y,z):
    fac = -g*cfg.M_mn / (x*x + y*y + (cfg.a_mn + np.sqrt(cfg.b_mn**2 + z*z))**2)**(1.5)
    return np.array([fac*x, fac*y, fac*z*(1.0+cfg.a_mn/np.sqrt(cfg.b_mn**2+z*z)) ])

 
# Exponential spheroid (bulge)
def mass_expon_spher(M_0,x):
	return M_0*(1.0-np.exp(-x)*(1.0+x+0.5*x**2))

def accel_expon_spher(a,rho_c,x,y,z):
    M_0 = 8.0*np.pi*a**3*rho_c
    r=np.sqrt(x*x+y*y+z*z)
    fac= -g*mass_expon_spher(M_0,r/a)/r**3
    return fac*np.array([x,y,z])


# Exponential thin disk (disk)
def dPhi_dR_efd(sigma_0,R_d,R,z):
    def integrand(R,z,k):
        return k*spec.jv(1,k*R)*np.exp(-k*np.abs(z))/(1.0+(k*R_d)**2)**(1.5)
    k_grid = np.linspace(0,10**2,10**5)
    integ = integrate.simps(integrand(R,z,k_grid),k_grid,even='avg')
    return 2.0*np.pi*g*sigma_0*R_d*R_d*integ

def dPhi_dz_efd(sigma_0,R_d,R,z):
    def integrand(R,z,k):
        return k*spec.jv(0,k*R)*np.exp(-k*np.abs(z))/(1.0+(k*R_d)**2)**(1.5)
    k_grid = np.linspace(0,10**2,10**5)
    integ = integrate.simps(integrand(R,z,k_grid),k_grid,even='avg')
    return 2.0*np.pi*g*sigma_0*R_d*R_d*integ*np.sign(z)  # sign(0)=0


def accel_expon_flat_disk(sigma_0,R_d,x,y,z):
    R=np.sqrt(x*x+y*y)
    aux=dPhi_dR_efd(sigma_0,R_d,R,z)
    if(R==0):
    	return -np.array([aux, aux, dPhi_dz_efd(sigma_0,R_d,R,z)])
    else:	
    	return -np.array([aux*(x/R), aux*(y/R), dPhi_dz_efd(sigma_0,R_d,R,z)])


# Exponential thin disk with spline interpolation 
def dPhi_dR_mesh(sigma_0,R_d,R_grid,z_grid):
    print('len(R_grid)=',len(R_grid))
    print('len(z_grid)=',len(z_grid))
    mesh = np.zeros((len(R_grid),len(z_grid)))
    for i in range(0,len(R_grid)):
        print('dPhi_dR: i=',i)
        for j in range(0,len(z_grid)):
            mesh[i,j]=dPhi_dR_efd(sigma_0,R_d,R_grid[i],z_grid[j])
    return mesh

def dPhi_dz_mesh(sigma_0,R_d,R_grid,z_grid):
    mesh = np.zeros((len(R_grid),len(z_grid)))
    for i in range(0,len(R_grid)):
        print('dPhi_dz: i=',i)
        for j in range(0,len(z_grid)):
            mesh[i,j]=dPhi_dz_efd(sigma_0,R_d,R_grid[i],z_grid[j])
    return mesh

R_grid = np.logspace(-10,2,1000)
z_grid = np.logspace(-10,1.7,1000)

# Uncomment the first time to evaluate
#mesh_R = dPhi_dR_mesh(cfg.sigma_efd,cfg.a_efd,R_grid,z_grid)
#mesh_z = dPhi_dz_mesh(cfg.sigma_efd,cfg.a_efd,R_grid,z_grid)

# Save
#np.savetxt('mesh_R.txt', mesh_R)
#np.savetxt('mesh_z.txt', mesh_z)

# Read
mesh_R = np.loadtxt('mesh_R.txt')
mesh_z = np.loadtxt('mesh_z.txt')

dPhi_dR_spline = interpolate.RectBivariateSpline(R_grid, z_grid, mesh_R,kx=3,ky=3)
dPhi_dz_spline = interpolate.RectBivariateSpline(R_grid, z_grid, mesh_z,kx=3,ky=3)


def accel_expon_flat_disk_spline(x,y,z):
    R=np.sqrt(x*x+y*y)
    abs_z =np.abs(z)
    return -np.array([dPhi_dR_spline(R,abs_z)[0][0]*x/R, dPhi_dR_spline(R,abs_z)[0][0]*y/R, np.sign(z)*dPhi_dz_spline(R,abs_z)[0][0]])


