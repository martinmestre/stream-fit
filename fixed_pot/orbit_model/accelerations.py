
# Gravitational accelerations
#-----------------------------



# Parameters
g= 4.300923924e-6


# Milky Way
#----------

# Dark matter

# RAR 
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
    fac= g*mass_rar_mw(r)/r**3
    return fac*np.array([x,y,z])


# NFW
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
    fac= g*mass_nfw_mw(r)/r**3
    return fac*np.array([x,y,z])


# Logarithmic
v_c = 220.0
q_phi = 0.9

def accel_logarithmic(x,y,z):
    r2 = x*x+y*y+z*z/(q_phi*q_phi)
    return v_c*v_c/r2*np.array([x,y,z/(q_phi*q_phi)])

# Barionic

M_gal = 2.325e7   # Unit of mass (in solar masses) used by Irrgang for G=1.
                  # Note that there G is not unity.

# Bulge (Plummer sphere)
M_b= 0.5e10
r_b= 0.23 
def accel_plum_bulge(x,y,z):
    fac = g*M_b/(r_b*r_b +x*x+y*y+z*z )**(1.5)
    return np.array([fac*x, fac*y, fac*z])

# Disk (Miyamoto-Nagai)
M_d=6.8e10  # Msun
a=3.0     # kpc
b=0.28    # kpc

def accel_mn_disk(x,y,z):
    fac = g*M_d / (x*x + y*y + (a + np.sqrt(b*b + z*z))**2)**(1.5)
    return np.array([fac*x, fac*y, fac*z*(1.0+a/np.sqrt(b*b+z*z)) ])

 