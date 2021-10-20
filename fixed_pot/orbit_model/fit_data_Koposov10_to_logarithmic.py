# Computation of the orbit model for Koposov+10 logarithmic potential and data.
# Optimization to find best fit orbit.


from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from scipy.integrate import solve_ivp
from astropy.io import ascii
from scipy import optimize
from scipy import interpolate
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import config as cfg
import accelerations_combine as acco
import GD1Koposov10_class as GD1_class


# Load Koposov data
exec(open("./stream_Koposov_data.py").read()) 

# Symplectic gradient
def symp_grad_mw(t, w):
    x=w[0]
    y=w[1]
    z=w[2]
    px=w[3]
    py=w[4]
    pz=w[5]
    r=np.sqrt(x*x+y*y+z*z)
    return [px,py,pz, acco.accel_mw(x,y,z)[0], acco.accel_mw(x,y,z)[1], acco.accel_mw(x,y,z)[2]]


# The orbit_model here defined works with cartesian coordinates at the input and sky-cartesian at output.


# Orbit model definition
#-------------------------------------
def orbit_model(w_0):
    print('w_0 = ',w_0)
    
    # ODE integration
    unit_t = 0.977792221680356   # Gyr
    time_span_s2 = 0.2/unit_t #
    t_0=0.0/unit_t    
    n_steps = 1000
    t_back = np.linspace(t_0,-time_span_s2, n_steps+1)
    t_forw = np.linspace(t_0,time_span_s2, n_steps+1)       
    sol_back = solve_ivp(symp_grad_mw, [t_0,-time_span_s2], w_0, t_eval=t_back,method='DOP853',rtol=5.0e-14,atol=0.5e-14)
    sol_forw = solve_ivp(symp_grad_mw, [t_0,time_span_s2], w_0, t_eval=t_forw,method='DOP853',rtol=5.0e-14,atol=0.5e-14)
    
    t = np.concatenate([sol_back.t,sol_forw.t])
    y = np.concatenate([sol_back.y, sol_forw.y],axis=1)
    y = np.delete(y,0,axis=1) #Remove duplicated column

    
    #Transformation to GD-1 frame of coordinates (\phi_1, \phi_2)
    galcen_distance = 8.5*u.kpc
    v_sun = coord.CartesianDifferential([11.1, 220.0+12.24, 7.25]*u.km/u.s)
    z_sun=0.0*u.kpc
    galac_coord=coord.Galactocentric(x=y[0]*u.kpc,y=y[1]*u.kpc,z=y[2]*u.kpc,
                                     v_x=y[3]*u.km/u.s,v_y=y[4]*u.km/u.s,v_z=y[5]*u.km/u.s,
                           galcen_distance=galcen_distance,galcen_v_sun=v_sun,z_sun=z_sun) 
    gd1_coord = galac_coord.transform_to(GD1_class.GD1Koposov10)
    phi_1 = gd1_coord.phi1
    phi_2 = gd1_coord.phi2
    d_hel = gd1_coord.distance
    v_hel = gd1_coord.radial_velocity
    mu_phi_1 = gd1_coord.pm_phi1_cosphi2/np.cos(phi_2)
    mu_phi_2 = gd1_coord.pm_phi2

    print(gd1_coord)
    return phi_1, phi_2, d_hel, v_hel, mu_phi_1, mu_phi_2,y[0],y[1],y[2] #the time "t" is not needed

# Test call:
#------------------
#w_0=np.array([-3.41, 13.0, 9.58, -200.4, -162.6, 13.9])
#w_0=np.array([-1.716384450138328077e+00, 1.401126977875537705e+01, 9.227965941367614278e+00,
#              -2.031152426286328421e+02, -1.368983530061972260e+02, 3.428826606652827280e+01])
w_0=np.array([-2.996988867007201574e+00, 1.309889922825600728e+01, 9.408043877712211511e+00,
              -2.006738213850691466e+02,-1.509052782505346215e+02, 2.239535561062263369e+01]) #Last fit

phi_1,phi_2,d_hel,v_hel,mu_1,mu_2,x,y,z = orbit_model(w_0) 
    

# Plot in galactocentric coordinates
#------------------------------------
fig=plt.figure(figsize=(10,10))

plt.scatter(x,y,s=0.1,marker='o',color='red')
plt.scatter(x,z,s=0.1,marker='o',color='blue')
plt.xlim(-15,10)
plt.ylim(-15,20)
plt.grid()
#plt.legend(fontsize='huge', handlelength=0.3)
#plt.tight_layout()
plt.show()
fig.savefig("plots/orbit_xz_fit_my_w0.png")

# Plots in the sky using the GD-1 frame
#------------------------------------------------------


# sky position
#fig=plt.figure(figsize=(10,7))
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4, 1, sharex=True, figsize=(7,20))


# Sky position
ax1.set_title('Koposov+10 data with Logarithmic potential')
ax1.scatter(phi_1.wrap_at(180*u.deg),phi_2,s=0.1,marker='o', color='red')
ax1.scatter(phi_1,phi_2,s=0.1,marker='o', color='violet')
ax1.errorbar(kop_sky['phi1'], kop_sky['phi2'], yerr=kop_sky['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')
ax1.set_ylim(-4,2)
ax1.set_ylabel(r'$\phi_2$ [degrees]')


# heliocentric radial velocity 
ax2.scatter(phi_1.wrap_at(180*u.deg),v_hel,s=0.1,marker='o', color='red')
ax2.errorbar(kop_rv['phi1'], kop_rv['vr'], yerr=kop_rv['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')
ax2.set_ylim(-400,200)
ax2.set_ylabel(r'$v_{\rm{LOS}}$ [km s$^{-1}$]')

# heliocentric distance
ax3.scatter(phi_1,d_hel,s=0.1,marker='o', color='violet')
ax3.errorbar(kop_dist['phi1'], kop_dist['dist'], yerr=kop_dist['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')
ax3.set_ylim(5,15)
ax3.set_ylabel(r'$d_\odot$ [kpc]')

# proper motion along phi_1 and phi_2
ax4.scatter(phi_1,mu_1,s=0.5, color='violet',label='$\mu_{\phi_1}$')
ax4.scatter(phi_1,mu_2,s=0.5, color='red',label='$\mu_{\phi_2}$')
ax4.errorbar(kop_pm['phi1'], kop_pm['mu_phi1'], yerr=kop_pm['err'], fmt='o', color='blue', label='Stream\n(Koposov+2010)')
ax4.errorbar(kop_pm['phi1'], kop_pm['mu_phi2'], yerr=kop_pm['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')
ax4.set_ylabel('$\mu$ [mas yr$^{-1}$]')
ax4.legend()

#plt.xlim(-80,20)
#plt.grid()
#plt.legend(fontsize='big', handlelength=0.3)
plt.xlabel(r'$\phi_1$ [degrees]')
plt.xlim(-150,100)
#plt.tight_layout()

plt.show()
fig.savefig("plots/sky_fit_kop10.png")

# Optimization
#-------------------------------------

def chi2(w_0):
    import wrap
    
    phi_1,phi_2,d_hel,v_hel,mu_1,mu_2,x,y,z = orbit_model(w_0) 
    cfg.phi_2_spl = interp1d(phi_1,phi_2,kind='cubic')
    cfg.d_hel_spl = interp1d(phi_1,d_hel,kind='cubic')
    cfg.v_hel_spl = interp1d(phi_1,v_hel,kind='cubic')
    cfg.mu_1_spl  = interp1d(phi_1,mu_1,kind='cubic')
    cfg.mu_2_spl  = interp1d(phi_1,mu_2,kind='cubic')


    cfg.phi_1_min = np.amin(phi_1.value)
    cfg.phi_1_max = np.amax(phi_1.value)

    
    sum=np.zeros(5)

    y_mod = wrap.phi_2_wrap(kop_sky['phi1'])
    y_dat = kop_sky['phi2']
    sigma2 = kop_sky['err']**2
    sum[0] = np.sum( (y_dat-y_mod)**2 / sigma2 )

    y_mod = wrap.v_hel_wrap(kop_rv['phi1'])
    y_dat = kop_rv['vr']
    sigma2 = kop_rv['err']**2
    sum[1] = np.sum( (y_dat-y_mod)**2 / sigma2 )

    y_mod = wrap.d_hel_wrap(kop_dist['phi1'])
    y_dat = kop_dist['dist']
    sigma2 = kop_dist['err']**2
    sum[2] = np.sum( (y_dat-y_mod)**2 / sigma2 )

    y_mod = wrap.mu_1_wrap(kop_pm['phi1'])
    y_dat = kop_pm['mu_phi1']
    sigma2 = kop_pm['err']**2
    sum[3] = np.sum( (y_dat-y_mod)**2 / sigma2 )

    y_mod = wrap.mu_2_wrap(kop_pm['phi1'])
    y_dat = kop_pm['mu_phi2']
    sigma2 = kop_pm['err']**2
    sum[4] = np.sum( (y_dat-y_mod)**2 / sigma2 )
    
    print('chi^2 =',np.sum(sum))
    return np.sum(sum)

print(chi2(w_0))


dx = 10.0
dv = 40.0
w_0=np.array([-3.41, 13.0, 9.58, -200.4, -162.6, 13.9]) 
bounds=((w_0[0]-dx,w_0[0]+dx), (w_0[1]-dx,w_0[1]+dx), (w_0[2]-dx,w_0[2]+dx),
        (w_0[3]-dv,w_0[3]+dv), (w_0[4]-dv,w_0[4]+dv), (w_0[5]-dv,w_0[5]+dv))

#opt=optimize.differential_evolution(chi2, bounds,strategy='best1bin',maxiter=20,popsize=20,tol=5.0e-8,atol=0.5e-8,disp=True,polish=True,workers=-1)

#param_fitted = opt.x

#np.savetxt('param_fit_Kop10.txt', param_fitted, delimiter=',')  
