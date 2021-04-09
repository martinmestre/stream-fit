# Multiple component accelerations

# Milky Way solution from Argüelles et al (20...)
#-------------------------------------------------

#Barionic component from Sofue 2013

# Inner bulge
def accel_bulge_inner(x,y,z):
	return accel_expon_spher(0.0038,3.6e13,x,y,z)

# Main bulge
def accel_bulge_main(x,y,z):
	return accel_expon_spher(0.12,1.9e11,x,y,z)

# Disk                    
R_d=3.0*0.93
m_d=4.4e10*1.2
sigma_0=m_d/(2.0*np.pi*R_d*R_d)
def accel_disk(x,y,z):
	#return accel_expon_flat_disk_spline(x,y,z)
	return accel_expon_flat_disk(sigma_0,R_d,x,y,z)

# Dark matter halo

def accel_halo_dm(x,y,z):
	return accel_rar_mw(x,y,z)

# Total contribution

def accel_mw(x,y,z):
	return accel_bulge_inner(x,y,z)+ accel_bulge_main(x,y,z)+ accel_disk(x,y,z) + accel_halo_dm(x,y,z)
	#return accel_logarithmic(x,y,z)
	#return accel_halo_dm(x,y,z)
	#return accel_disk(x,y,z)
