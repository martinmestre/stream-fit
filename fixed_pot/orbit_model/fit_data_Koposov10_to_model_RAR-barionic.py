# Computation of the orbit model for RAR + barionic component,
# using Koposov+20 data.
# Optimization to find best fit orbit.

exec(open("./stream_Koposov_data.py").read()) 

def orbit_model(alpha,delta,distance,mu_alpha,mu_delta,v_los):
    print('param= ',alpha,delta,distance,mu_alpha,mu_delta,v_los)
    
    # Transformation to galactocentric coordinates
    sky_coord = coord.ICRS(ra=alpha*u.degree, dec=delta*u.degree,
                distance=distance*u.kpc,
                pm_ra_cosdec=mu_alpha*np.cos(delta*u.degree)*u.mas/u.yr,
                pm_dec=mu_delta*u.mas/u.yr,
                radial_velocity=v_los*u.km/u.s)
    galcen_distance = 8.129*u.kpc
    v_sun = coord.CartesianDifferential([11.1, 229.0+12.24, 7.25]*u.km/u.s)
    z_sun=0.0*u.kpc
    frame = coord.Galactocentric(galcen_distance=galcen_distance,
                                galcen_v_sun=v_sun,
                                z_sun=z_sun)
    galac_coord= sky_coord.transform_to(frame)
    
    w_0 = np.zeros(6)
    w_0[:3]=[galac_coord.x/u.kpc,galac_coord.y/u.kpc,galac_coord.z/u.kpc]
    w_0[3:]=[galac_coord.v_x/(u.km/u.s),galac_coord.v_y/(u.km/u.s),galac_coord.v_z/(u.km/u.s)]

    
    # ODE integration
    unit_t = 0.977792221680356   # Gyr
    time_span_s2 = 0.15/unit_t #
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
    galac_coord=coord.Galactocentric(x=y[0]*u.kpc,y=y[1]*u.kpc,z=y[2]*u.kpc,
                                     v_x=y[3]*u.km/u.s,v_y=y[4]*u.km/u.s,v_z=y[5]*u.km/u.s,
                           galcen_distance=galcen_distance,galcen_v_sun=v_sun,z_sun=z_sun) 
    gd1_coord = galac_coord.transform_to(GD1Koposov10)
    phi_1 = gd1_coord.phi1
    phi_2 = gd1_coord.phi2
    d_hel = gd1_coord.distance
    v_hel = gd1_coord.radial_velocity
    mu_phi_1 = gd1_coord.pm_phi1_cosphi2/np.cos(phi_2)  #not used by Ibata
    mu_phi_2 = gd1_coord.pm_phi2
    return phi_1, phi_2, d_hel, v_hel, mu_phi_1, mu_phi_2, y[0], y[1], y[2]
    # Transformation to ICRS coordinates      
    #icrs_coord=galac_coord.transform_to(coord.ICRS)
    #mu_ra = icrs_coord.pm_ra_cosdec / np.cos(icrs_coord.dec)
    #mu_dec= icrs_coord.pm_dec
    #return t, phi_1, phi_2, d_hel, v_hel, mu_ra, mu_dec

# test call:
valor_pm_ra=-6.53  #adimensional
valor_dec=43.717
valor_cosdec= np.cos(valor_dec*u.deg)
valor_pm_ra_cosdec=valor_pm_ra*valor_cosdec
# This is the initial value I used for optimizer:
#phi_1,phi_2,d_hel,v_hel,mu_1,mu_2 = orbit_model(157.6,valor_dec, 8.25, valor_pm_ra_cosdec ,-11.0 , -90.0) 
# This is the initial value I used the second time I refined the search:
#phi_1,phi_2,d_hel,v_hel,mu_1,mu_2 = orbit_model(154.43815737978053, 41.352611518783334, 9.127975838807465, -6.135527680412597, -9.275558167001593, -63.00000001)  
# This are the fitter orbit parameters:
phi_1,phi_2,d_hel,v_hel,mu_1,mu_2,x,y,z = orbit_model(1.537274943846678070e+02, 4.059754233752208563e+01, 7.854749086400735436e+00, -7.976185984536376061e+00, -1.168788612180591713e+01, -4.410000000700000555e+01)

print(phi_1)
print(d_hel)

# Plot in galactocentric coordinates
#------------------------------------
fig=plt.figure(figsize=(10,10))

plt.scatter(x,y,s=0.1,marker='o',color='red')
plt.scatter(x,z,s=0.1,marker='o',color='blue')
plt.xlim(-15,10)
plt.ylim(-15,20)

# Plots in the sky using the GD-1 frame
#---------------------------------------

# sky position
#plt.scatter(phi_1.wrap_at(180*u.deg),phi_2,s=0.1,marker='o', color='red')
#plt.scatter(phi_1,phi_2,s=0.1,marker='o', color='violet')
#plt.errorbar(kop_sky['phi1'], kop_sky['phi2'], yerr=kop_sky['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')
#plt.ylim(-4,2)

# heliocentric radial velocity 
#plt.scatter(phi_1.wrap_at(180*u.deg),v_hel,s=0.1,marker='o', color='red')
#plt.errorbar(kop_rv['phi1'], kop_rv['vr'], yerr=kop_rv['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')
#plt.ylim(-400,200)

# heliocentric distance
#plt.scatter(phi_1,d_hel,s=0.1,marker='o', color='violet')
#plt.errorbar(kop_dist['phi1'], kop_dist['dist'], yerr=kop_dist['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')
#plt.ylim(5,15)

# proper motion along phi_1 and phi_2
#plt.scatter(phi_1,mu_1,s=0.1,marker='o', color='violet')
#plt.scatter(phi_1,mu_2,s=0.1,marker='o', color='red')
#plt.errorbar(kop_pm['phi1'], kop_pm['mu_phi1'], yerr=kop_pm['err'], fmt='o', color='blue', label='Stream\n(Koposov+2010)')
#plt.errorbar(kop_pm['phi1'], kop_pm['mu_phi2'], yerr=kop_pm['err'], fmt='o', color='cyan', label='Stream\n(Koposov+2010)')

#plt.xlim(-80,20)
plt.grid()
plt.legend(fontsize='small', handlelength=0.3)


plt.tight_layout()
plt.show()
fig.savefig("plots/orbit_xz_fit_data_Kop2rar.png")


# Optimization
#--------------------

def chi2(w_0):
    phi_1,phi_2,d_hel,v_hel,mu_1,mu_2,x,y,z = orbit_model(w_0[0],w_0[1],w_0[2],w_0[3],w_0[4],w_0[5]) 
    phi_2_spl = interp1d(phi_1,phi_2,kind='cubic')
    d_hel_spl = interp1d(phi_1,d_hel,kind='cubic')
    v_hel_spl = interp1d(phi_1,v_hel,kind='cubic')
    mu_1_spl  = interp1d(phi_1,mu_1,kind='cubic')
    mu_2_spl  = interp1d(phi_1,mu_2,kind='cubic')

    phi_1_min = np.amin(phi_1.value)
    phi_1_max = np.amax(phi_1.value)
    #Wrappers
    def phi_2_wrap(x):
        array = np.zeros(len(x))
        for i in range(0,len(x)):
            if(x[i]<phi_1_min):
                array[i]= phi_2_spl(phi_1_min)
            elif(x[i]>phi_1_max):
                array[i]= phi_2_spl(phi_1_max)
            else:
                array[i]= phi_2_spl(x[i])
        return array

    def d_hel_wrap(x):
        array = np.zeros(len(x))
        for i in range(len(x)):
            if(x[i]<phi_1_min):
                array[i]= d_hel_spl(phi_1_min)
            elif(x[i]>phi_1_max):
                array[i]= d_hel_spl(phi_1_max)
            else:
                array[i]= d_hel_spl(x[i])
        return array

    def v_hel_wrap(x):
        array = np.zeros(len(x))
        for i in range(0,len(x)):
            if(x[i]<phi_1_min):
                array[i]= v_hel_spl(phi_1_min)
            elif(x[i]>phi_1_max):
                array[i]= v_hel_spl(phi_1_max)
            else:
                array[i]= v_hel_spl(x[i])
        return array

    def mu_1_wrap(x):
        array = np.zeros(len(x))
        for i in range(0,len(x)):
            if(x[i]<phi_1_min):
                array[i]= mu_1_spl(phi_1_min)
            elif(x[i]>phi_1_max):
                array[i]= mu_1_spl(phi_1_max)
            else:
                array[i]= mu_1_spl(x[i])
        return array
    
    def mu_2_wrap(x):
        array = np.zeros(len(x))
        for i in range(0,len(x)):
            if(x[i]<phi_1_min):
                array[i]= mu_2_spl(phi_1_min)
            elif(x[i]>phi_1_max):
                array[i]= mu_2_spl(phi_1_max)
            else:
                array[i]= mu_2_spl(x[i])
        return array
        
    sum=np.zeros(5)
    y_mod = phi_2_wrap(kop_sky['phi1'])
    y_dat = kop_sky['phi2']
    sigma2 = kop_sky['err']**2
    sum[0] = np.sum( (y_dat-y_mod)**2 / sigma2 )

    y_mod = v_hel_wrap(kop_rv['phi1'])
    y_dat = kop_rv['vr']
    sigma2 = kop_rv['err']**2
    sum[1] = np.sum( (y_dat-y_mod)**2 / sigma2 )
    
    y_mod = d_hel_wrap(kop_dist['phi1'])
    y_dat = kop_dist['dist']
    sigma2 = kop_dist['err']**2
    sum[2] = np.sum( (y_dat-y_mod)**2 / sigma2 )
    
    y_mod = mu_1_wrap(kop_pm['phi1'])
    y_dat = kop_pm['mu_phi1']
    sigma2 = kop_pm['err']**2
    sum[3] = np.sum( (y_dat-y_mod)**2 / sigma2 )
    
    y_mod = mu_2_wrap(kop_pm['phi1'])
    y_dat = kop_pm['mu_phi2']
    sigma2 = kop_pm['err']**2
    sum[4] = np.sum( (y_dat-y_mod)**2 / sigma2 )
    
    print('chi^2 =',np.sum(sum))
    return np.sum(sum)




valor_pm_ra=-6.53  #adimensional
valor_dec=43.717
valor_cosdec= np.cos(valor_dec*u.deg)
#w_0=np.array([154.43815737978053, 41.352611518783334, 9.127975838807465, -6.135527680412597, -9.275558167001593, -63.00000001]) 
w_0= np.array([1.537274943846678070e+02, 4.059754233752208563e+01, 7.854749086400735436e+00, -7.976185984536376061e+00, -1.168788612180591713e+01, -4.410000000700000555e+01])
dw = np.abs(w_0)*0.3
bounds=((w_0[0]-dw[0],w_0[0]+dw[0]), (w_0[1]-dw[1],w_0[1]+dw[1]), (w_0[2]-dw[2],w_0[2]+dw[2]),
        (w_0[3]-dw[3],w_0[3]+dw[3]), (w_0[4]-dw[4],w_0[4]+dw[4]), (w_0[5]-dw[5],w_0[5]+dw[5]))

#opt=optimize.differential_evolution(chi2, bounds,strategy='best1bin',maxiter=20,popsize=20,tol=5.0e-8,atol=0.5e-8,disp=True,polish=True,workers=-1)

#param_fitted = opt.x

#np.savetxt('param_fitted.txt', param_fitted, delimiter=',')  

print('chi2(w_0)=',chi2(w_0))