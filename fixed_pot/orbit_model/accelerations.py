
# Gravitational accelerations
#-----------------------------
import scipy.special as spec


# Parameters
g= 4.300923924e-6


# Milky Way
#----------



# RAR (halo)
def rho_rar_mw(r):
    if(r<r_max_mw):
        return rho_spline(r)
    else:
        return rho_spline(r_max_mw)

def mass_rar_mw(r):
    if(r<r_max_mw):
        return mass_spline(r)
    else:
        return mass_spline(r_max_mw)
    
def accel_rar_mw(x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    fac= -g*mass_rar_mw(r)/r**3
    return fac*np.array([x,y,z])


# NFW (halo)
H_0 = 0.0674
rho_c =3.0*H_0**2/(8.0*np.pi*g)
c = 12.
m_200 = 0.97e12
r_200 = (3.0*m_200/(800.0*np.pi*rho_c))**(1.0/3.0)
a_nfw = r_200/c
print('r_200=',r_200)
print('a_nfw=',a_nfw)

def f(x):
    return np.log(1.0+x)- x/(1.0+x)

def rho_nfw_mw(r):
    rho_0 = m_200/(4.0*np.pi*f(c)*a_nfw**3)
    return rho_0/(r/a_nfw)/(1+r/a_nfw)**2

def mass_nfw_mw(r):
    return (m_200/f(c))*f(r/a_nfw)


def accel_nfw_mw(x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    fac= -g*mass_nfw_mw(r)/r**3
    return fac*np.array([x,y,z])


# Logarithmic (halo)
v_c = 220.0
q_phi = 0.9

def accel_logarithmic(x,y,z):
    r2 = x*x+y*y+z*z/(q_phi*q_phi)
    return -v_c*v_c/r2*np.array([x,y,z/(q_phi*q_phi)])




# Plummer sphere (bulge)
M_b= 0.5e10
r_b= 0.23 
def accel_plum_bulge(x,y,z):
    fac = -g*M_b/(r_b*r_b +x*x+y*y+z*z )**(1.5)
    return np.array([fac*x, fac*y, fac*z])

# Miyamoto-Nagai (disk)
M_d=6.8e10  # Msun
a=3.0     # kpc
b=0.28    # kpc

def accel_mn_disk(x,y,z):
    fac = -g*M_d / (x*x + y*y + (a + np.sqrt(b*b + z*z))**2)**(1.5)
    return np.array([fac*x, fac*y, fac*z*(1.0+a/np.sqrt(b*b+z*z)) ])

 
# Exponential spheroid (bulge)
def mass_expon_spher(a,rho_c,x):
	M_0 = 8.0*np.pi*a**3*rho_c
	return M_0*(1.0-np.exp(-x)*(1.0+x+0.5*x**2))

def accel_expon_spher(a,rho_c,x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    fac= -g*mass_expon_spher(a,rho_c,r/a)/r**3
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
    if(z==0):
    	return 0.0
    else:
    	return 2.0*np.pi*g*sigma_0*R_d*R_d*integ*np.sign(z)  


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
        for j in range(0,len(z_grid)):
            mesh[i,j]=dPhi_dR_efd(sigma_0,R_d,R_grid[i],z_grid[j])
    return mesh

def dPhi_dz_mesh(sigma_0,R_d,R_grid,z_grid):
    mesh = np.zeros((len(R_grid),len(z_grid)))
    for i in range(0,len(R_grid)):
        for j in range(0,len(z_grid)):
            mesh[i,j]=dPhi_dz_efd(sigma_0,R_d,R_grid[i],z_grid[j])
    return mesh

R_d=3.0
m_d=4.4e10
sigma_0=m_d/(2.0*np.pi*R_d*R_d)
R_grid = np.linspace(0,30.0,200)
z_grid = np.linspace(-15.0,15.0,200)

# Uncomment the first time to evaluate
#mesh_R = dPhi_dR_mesh(sigma_0,R_d,R_grid,z_grid)
#mesh_z = dPhi_dz_mesh(sigma_0,R_d,R_grid,z_grid)

# Read
mesh_R = np.loadtxt('mesh_R.txt')
mesh_z = np.loadtxt('mesh_z.txt')

# Save
#np.savetxt('mesh_R.txt', mesh_R)
#np.savetxt('mesh_z.txt', mesh_z)

dPhi_dR_spline = interpolate.RectBivariateSpline(R_grid, z_grid, mesh_R,kx=3,ky=3)
dPhi_dz_spline = interpolate.RectBivariateSpline(R_grid, z_grid, mesh_z,kx=3,ky=3)

def accel_expon_flat_disk_spline(x,y,z):
	R=np.sqrt(x*x+y*y)
	return -np.array([dPhi_dR_spline(R,z)[0][0]*x/R, dPhi_dR_spline(R,z)[0][0]*y/R, dPhi_dz_spline(R,z)[0][0]])

#print ('accel_efds=',accel_expon_flat_disk_spline(1.,1.,1.))
#print('spline=',dPhi_dR_spline(np.sqrt(2.),1.)[0])
#print('spline=',dPhi_dR_spline(np.sqrt(2.),1.)[0][0])
#print('fun=',dPhi_dR_efd(sigma_0,R_d,np.sqrt(2.),1.))
