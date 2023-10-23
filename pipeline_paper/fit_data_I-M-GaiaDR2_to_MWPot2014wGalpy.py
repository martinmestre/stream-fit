
# Computation of the orbit model the best fit potential in Malhan+19 (MWPotential2014+axisymmetry),
# using Malhan+19 data.
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
import GD1Koposov10_class as GD1_class
import pandas as pd
from astroquery.gaia import Gaia
import galpy as gp
from galpy.potential.mwpotentials import MWPotential2014
from galpy.potential import TriaxialNFWPotential,PowerSphericalPotentialwCutoff,MiyamotoNagaiPotential
from galpy.orbit import Orbit
from galpy.util import conversion

# Orbit definition
# The orbit_model here defined works with sky coordinates at input and sky-cartesian at output
def orbit_model(alpha,delta,distance,mu_alpha_cosdelta,mu_delta,v_los):
#    print('param= ',alpha,delta,distance,mu_alpha_cosdelta,mu_delta,v_los)

    # Transformation to galactocentric coordinates
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import CartesianDifferential
    c= SkyCoord(ra=alpha*u.degree, dec=delta*u.degree,
                distance=distance*u.kpc,
                pm_ra_cosdec=mu_alpha_cosdelta*u.mas/u.yr,
                pm_dec=mu_delta*u.mas/u.yr,
                radial_velocity=v_los*u.km/u.s,
				galcen_distance = 8.0*u.kpc,
				galcen_v_sun = coord.CartesianDifferential([11.1, v_circ_sun+12.24, 7.25]*u.km/u.s),
				z_sun=0.0*u.kpc)

    o=Orbit(c)
    n_steps=1000
    # Backward integration
    ts1= np.linspace(0,-0.2,n_steps+1)*u.Gyr
    o.integrate(ts1,mw,method='dop853_c')
    sol_back = np.array([o.ra(ts1),o.dec(ts1),o.dist(ts1),o.pmra(ts1),o.pmdec(ts1),o.vlos(ts1)])
    ts2= np.linspace(0,0.2,n_steps+1)*u.Gyr
    o.integrate(ts2,mw,method='dop853_c')
    sol_forw = np.array([o.ra(ts2),o.dec(ts2),o.dist(ts2),o.pmra(ts2),o.pmdec(ts2),o.vlos(ts2)])

    t = np.concatenate([ts1,ts2])
    y_back = np.delete(sol_back, 0, axis=1)  # Remove duplicated column
    y_back = np.flip(y_back, axis=1)
    y = np.concatenate([y_back, sol_forw],  axis=1)
    y = np.flip(y, axis=1)


    #Transformation to GD-1 frame of coordinates (\phi_1, \phi_2)
    icrs_coord=coord.ICRS(ra=y[0]*u.degree, dec=y[1]*u.degree,
                distance=y[2]*u.kpc,
                pm_ra_cosdec=y[3]*u.mas/u.yr,
                pm_dec=y[4]*u.mas/u.yr,
                radial_velocity=y[5]*u.km/u.s)
    mu_ra_cosdec = icrs_coord.pm_ra_cosdec
    mu_dec= icrs_coord.pm_dec
    gd1_coord = icrs_coord.transform_to(GD1_class.GD1Koposov10)
    phi_1 = gd1_coord.phi1
    phi_2 = gd1_coord.phi2
    d_hel = gd1_coord.distance
    v_hel = gd1_coord.radial_velocity
    #mu_phi_1 = gd1_coord.pm_phi1_cosphi2/np.cos(phi_2)  #not used by Ibata
    #mu_phi_2 = gd1_coord.pm_phi2
    #return phi_1, phi_2, d_hel, v_hel, mu_phi_1, mu_phi_2

    np.savetxt('observable_orbit_NFW-MW.txt', (phi_1, phi_2, d_hel, mu_ra_cosdec, mu_dec, v_hel))
    return phi_1, phi_2, d_hel, mu_ra_cosdec, mu_dec, v_hel


# Ibata's polynomials.
# [x] = radians
# [S] = radians
# [D] = kpc
# [V] = km/s
# [MU] = mas/year
class IbaPoly:
    def S(self,x):
        return 0.008367*x**3-0.05332*x**2-0.07739*x-0.02007
    def D(self,x):
        return -4.302*x**5-11.54*x**4-7.161*x**3+5.985*x**2+8.595*x+10.36
    def V(self,x):
        return 90.68*x**3+204.5*x**2-254.2*x-261.5
    def MU_RA(self,x):
        return 3.794*x**3+9.467*x**2+1.615*x-7.844
    def MU_DEC(self,x):
        return -1.225*x**3+8.313*x**2+18.68*x-3.95
    limit = [-90,10]

pol = IbaPoly()
print('valor=',pol.S(np.linspace(-90,10,10)))

# Observations (Ibata polynomials evaluated in a grid)
Iba_sky = pd.DataFrame()
Iba_sky['phi_1'] = np.linspace(IbaPoly.limit[0],IbaPoly.limit[1],100)
Iba_sky['phi_2'] = pol.S(cfg.deg2rad*Iba_sky['phi_1'])/cfg.deg2rad
Iba_sky['d_hel'] = pol.D(cfg.deg2rad*Iba_sky['phi_1'])
Iba_sky['v_hel'] = pol.V(cfg.deg2rad*Iba_sky['phi_1'])
Iba_sky['mu_ra'] = pol.MU_RA(cfg.deg2rad*Iba_sky['phi_1'])
Iba_sky['mu_dec'] = pol.MU_DEC(cfg.deg2rad*Iba_sky['phi_1'])

def invert_ic(u_0):
    w_0=u_0
    gd1_coord = GD1_class.GD1Koposov10(phi1=u_0[0]*u.degree, phi2=u_0[1]*u.degree)
    icrs_coord = gd1_coord.transform_to(coord.ICRS)
    w_0[0]=icrs_coord.ra.value
    w_0[1]=icrs_coord.dec.value
    return w_0

#print('To compare with the Figures in Malhan+19')
# [Iba_sky['RA'],Iba_sky['Dec']] = invert_ic([Iba_sky['phi_1'],Iba_sky['phi_2']])


# Load Malhan+19 tables
"""
Uncomment to download data from Gaia for the first time.
Otherwise just open the file at the bottom of this cell.

data = cross_match('asu.tsv')
"""

# Open cross matched file
# from astropy.table import Table
# data = Table.read('cross_matched.dat', format='ascii')
# print(type(data))
# print(len(data))

# arr = data['Pmember']
# boolarr=  (arr=='Y')

# data = data[boolarr]
# print(type(data))
# print(len(data))



# Chi^2 function
def chi2(w_0):
	import wrap_ra as wrap

	phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel = orbit_model(w_0[0],w_0[1],w_0[2],w_0[3],w_0[4],w_0[5])

	[ra,dec] = invert_ic([phi_1.value,phi_2.value])
	cfg.dec_spl = interp1d(ra,dec,kind='cubic')
	cfg.v_hel_spl = interp1d(ra,v_hel,kind='cubic')
	cfg.mu_ra_spl  = interp1d(ra,mu_ra,kind='cubic')
	cfg.mu_dec_spl  = interp1d(ra,mu_dec,kind='cubic')

	cfg.ra_min = np.amin(ra)
	cfg.ra_max = np.amax(ra)

	sum=np.zeros(4)
	y_mod = wrap.dec_wrap(data['RAJ2000'])
	y_dat = data['DEJ2000']
	sigma2 = 0.1**2
	sum[0] = np.sum( (y_dat-y_mod)**2 / sigma2 )

	y_mod = wrap.mu_ra_wrap(data['RAJ2000'])
	y_dat = data['pmra']
	sigma2 = data['pmra_error']**2
	sum[1] = np.sum( (y_dat-y_mod)**2 / sigma2 )

	y_mod = wrap.mu_dec_wrap(data['RAJ2000'])
	y_dat = data['pmdec']
	sigma2 = data['pmdec_error']**2
	sum[2] = np.sum( (y_dat-y_mod)**2 / sigma2 )

	y_mod = wrap.v_hel_wrap(data['RAJ2000'])
	y_dat = data['Vlos']
	sigma2 = data['e_Vlos']**2
	sum[3] = np.sum( (y_dat-y_mod)**2 / sigma2 )

	print('chi^2 =',np.sum(sum))
	return np.sum(sum)

sigma_array = np.array([0.5, 1.5, 2.0, 2.0, 10.0])
def chi2_stream(w_0):
    """Chi^2 stream function."""
    import wrap

    phi_1, phi_2, d_hel, mu_ra, mu_dec, v_hel = orbit_model(w_0[0],w_0[1],w_0[2],w_0[3],w_0[4],w_0[5])
    cfg.phi_2_spl = interp1d(phi_1, phi_2, kind='cubic')
    cfg.d_hel_spl = interp1d(phi_1, d_hel, kind='cubic')
    cfg.v_hel_spl = interp1d(phi_1, v_hel,  kind='cubic')
    cfg.mu_ra_spl = interp1d(phi_1, mu_ra, kind='cubic')
    cfg.mu_dec_spl = interp1d(phi_1, mu_dec, kind='cubic')
    cfg.phi_1_min = np.amin(phi_1.value)
    cfg.phi_1_max = np.amax(phi_1.value)

    sum = np.zeros(5)

    y_mod = wrap.phi_2_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['phi_2']
    sigma2 = sigma_array[0]**2
    sum[0] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.d_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['d_hel']
    sigma2 = sigma_array[1]**2
    sum[1] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_ra_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_ra']
    sigma2 = sigma_array[2]**2
    sum[2] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.mu_dec_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['mu_dec']
    sigma2 = sigma_array[3]**2
    sum[3] = np.sum((y_dat-y_mod)**2 / sigma2)

    y_mod = wrap.v_hel_wrap(Iba_sky['phi_1'])
    y_dat = Iba_sky['v_hel']
    sigma2 = sigma_array[4]**2
    sum[4] = np.sum((y_dat-y_mod)**2 / sigma2)

    print('chi^2_stream =',np.sum(sum))
    return np.sum(sum)


# Optimization
#----------------

# Potential selection:

bp= PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
mp= MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
hp=TriaxialNFWPotential(a=16.0*u.kpc,b=1.0,c=0.82,normalize=0.59)
mw = bp+mp+hp
print('mw=',mw)
# for i in range(0,3):
#    mw[i].turn_physical_on()

# r_sun=8.0*u.kpc
r_sun=8.122*u.kpc
v0=mw[0].vcirc(r_sun)
v1=mw[1].vcirc(r_sun)
v2=mw[2].vcirc(r_sun)
vo=220.0
print("v_i=",v0,v1,v2)
v_circ_sun=np.sqrt(v0*v0+v1*v1+v2*v2)*vo
print('v_circ_sun=',v_circ_sun)
a_R = np.zeros(3)
a_z = np.zeros(3)
r = np.array([14.0, 3.0])/8.0
for i in range(0,3):
    a_R[i] = mw[i].Rforce(r[0],r[1])*conversion.force_in_kmsMyr(220.,8.)
    a_z[i] = mw[i].zforce(r[0],r[1])*conversion.force_in_kmsMyr(220.,8.)
    print("i a_R a_z = ", i,  a_R[i], a_z[i])

a = np.sqrt([np.dot(a_R,a_R), np.dot(a_z,a_z)])
a_mod = np.sqrt(np.dot(a,a))
print("a=",a, "  a_mod=",a_mod)

# Initial conditions:
k0=50
u_0=np.array([Iba_sky['phi_1'][k0],Iba_sky['phi_2'][k0],Iba_sky['d_hel'][k0],
              Iba_sky['mu_ra'][k0],Iba_sky['mu_dec'][k0],Iba_sky['v_hel'][k0]])
w_0 = invert_ic(u_0)
print('We used seed: w_0=',w_0)

dw = np.abs(w_0)*0.5
dw[0]=0.1
bounds=((w_0[0]-dw[0],w_0[0]+dw[0]), (w_0[1]-dw[1],w_0[1]+dw[1]), (w_0[2]-dw[2],w_0[2]+dw[2]),
        (w_0[3]-dw[3],w_0[3]+dw[3]), (w_0[4]-dw[4],w_0[4]+dw[4]), (w_0[5]-dw[5],w_0[5]+dw[5]))
# bounds =((130,220),(20,60),(7,12),(-20,0),(-14,-1),(-300,250) )
# bounds =((0,360),(-90,90),(4,12),(-20,0),(-20,0),(-70,-20) )
# A posteriori refinement:
# w_0 = np.array([1.493370985649168858e+02, 3.669966976308609219e+01, 7.917039545144660018e+00,
#                 -7.050282547954606294e+00, -1.254565799483599520e+01, -1.636083097847286538e+01])
# dw = np.abs(w_0)*1.e-1
# bounds = ((w_0[0]-dw[0], w_0[0]+dw[0]), (w_0[1]-dw[1], w_0[1]+dw[1]), (w_0[2]-dw[2], w_0[2]+dw[2]),
#           (w_0[3]-dw[3], w_0[3]+dw[3]), (w_0[4]-dw[4], w_0[4]+dw[4]), (w_0[5]-dw[5], w_0[5]+dw[5]))

# opt=optimize.differential_evolution(chi2_stream, bounds,strategy='best2bin',maxiter=200,popsize=200,
#                                      tol=5.0e-8,atol=0.0,disp=True,polish=True,workers=-1)

# param_fitted = opt.x

# np.savetxt('param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt', param_fitted, delimiter=',')


# # Test call:
w_0 = np.loadtxt('param_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.txt')



print("w_0=",w_0)
phi_1,phi_2,d_hel,mu_ra,mu_dec,v_hel = orbit_model(w_0[0],w_0[1],w_0[2],w_0[3],w_0[4],w_0[5])
[ra,dec] = invert_ic([phi_1.value,phi_2.value])
chi2_stream(w_0)



# Plots in the sky using the GD-1 frame
#---------------------------------------
fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex=True, figsize=(10,15))
font = {"size": 20}
plt.rc('font', **font)

# Sky position
# ax1.set_title('Axisym. NFW in MWPotential2014 fit by Malhan+19')
ax1.scatter(phi_1.wrap_at(180*u.deg), phi_2, s=0.1, marker='o', color='red',label='Fit (Malhan+19)' )
ax1.plot(Iba_sky['phi_1'], Iba_sky['phi_2'], color='blue', label='Data (Ibata+20)')
ax1.set_ylim(-4, 2)
ax1.set_ylabel(r'$\phi_2$ [degrees]',fontsize=20)
ax1.legend()

# Heliocentric radial velocity
ax2.scatter(phi_1.wrap_at(180*u.deg), v_hel, s=0.1, marker='o', color='red')
ax2.plot(Iba_sky['phi_1'], Iba_sky['v_hel'], color='blue')
ax2.set_ylim(-300, 300)
ax2.set_ylabel(r'$v_{\rm{LOS}}$ [km s$^{-1}$]',fontsize=20)

# Heliocentric distance
ax3.scatter(phi_1, d_hel, s=0.1, marker='o', color='red')
ax3.plot(Iba_sky['phi_1'], Iba_sky['d_hel'], color='blue')
ax3.set_ylim(7, 12)
ax3.set_ylabel(r'$D$ [kpc]', fontsize=20)

# Proper motion along RA
ax4.scatter(phi_1, mu_ra, s=0.5, color='red')
ax4.plot(Iba_sky['phi_1'], Iba_sky['mu_ra'], color='blue')
ax4.set_ylabel(r'$\mu_\alpha$ [mas yr$^{-1}$]', fontsize=20)
ax4.set_ylim(-10, 0)


# Proper motion along DEC
ax5.scatter(phi_1, mu_dec, s=0.5, color='red')
ax5.plot(Iba_sky['phi_1'], Iba_sky['mu_dec'], color='blue')
ax5.set_ylabel(r'$\mu_\delta$ [mas yr$^{-1}$]', fontsize=20)


plt.xlabel(r'$\phi_1$ [degrees]', fontsize=20)
plt.xlim(IbaPoly.limit[0], IbaPoly.limit[1])
plt.tight_layout()
fig.savefig("plots/sky_fit_I-M-GaiaDR2_to_MWPot2014wGalpy_Polys.png")


# Plots using RA
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, sharex=True, figsize=(13,15))

# proper motion along RA
ax1.scatter(ra,mu_ra,s=0.5,color='red')
ax1.scatter(ra,mu_dec,s=0.5,color='orange')
ax1.plot(Iba_sky['RA'], Iba_sky['mu_ra'], color='blue', label=r'$\mu_\alpha$(Malhan+19)')
ax1.plot(Iba_sky['RA'], Iba_sky['mu_dec'], color='orange', label=r'$\mu_\delta$(Malhan+19)')
ax1.errorbar(data['RAJ2000'],data['pmra'],yerr=data['pmra_error'],fmt='o',color='black')
ax1.errorbar(data['RAJ2000'],data['pmdec'],yerr=data['pmdec_error'],fmt='o',color='black')
ax1.set_ylabel(r'$\mu$ [mas yr$^{-1}$]')
ax1.set_ylim(-15,5)
ax1.legend()

# heliocentric radial velocity
ax2.scatter(ra,v_hel,s=0.5,color='red')
ax2.plot(Iba_sky['RA'], Iba_sky['v_hel'], color='blue', label='Ibata+19')
ax2.errorbar(data['RAJ2000'],data['Vlos'],yerr=data['e_Vlos'],fmt='o',color='black')
ax2.set_ylim(-300,300)
ax2.set_ylabel(r'$v_{\rm{LOS}}$ [km s$^{-1}$]')
ax2.set_xlim(130,230)

# sky position
ax3.scatter(ra,dec,s=0.5,color='red')
ax3.plot(Iba_sky['RA'], Iba_sky['Dec'], color='blue', label='Ibata+19')
ax3.scatter(data['RAJ2000'],data['DEJ2000'])
ax3.set_ylim(10,70)
ax3.set_ylabel(r'$\delta$ [degrees]')
ax3.set_xlim(130,230)

plt.xlabel('RA [degrees]')
plt.xlim(130,230)
plt.tight_layout()
plt.show()
fig.savefig("plots/sky_fit_I-M-GaiaDR2_to_MWPot2014wGalpy.png")




