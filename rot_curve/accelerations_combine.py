# Multiple component accelerations

# Milky Way solution from Argüelles et al (20...)
#-------------------------------------------------

import accelerations as ac
import numpy as np
import config as cfg

#Barionic component from Sofue 2013

# Inner bulge
def accel_bulge_inner(x,y,z):
	return ac.accel_expon_spher(cfg.a_exps_in,cfg.rho_exps_in,x,y,z)

# Main bulge
def accel_bulge_main(x,y,z):
	return ac.accel_expon_spher(cfg.a_exps_ma,cfg.rho_exps_ma,x,y,z)

# Disk                    

def accel_disk(x,y,z):
	return ac.accel_expon_flat_disk_spline(x,y,z)
	#return ac.accel_expon_flat_disk(cfg.sigma_efd,cfg.a_efd,x,y,z)

# Dark matter halo

def accel_halo_dm(x,y,z):
	return ac.accel_rar_mw(x,y,z)

# Total contribution

def accel_mw(x,y,z):
	return accel_bulge_inner(x,y,z)+ accel_bulge_main(x,y,z)+ accel_disk(x,y,z) + accel_halo_dm(x,y,z)
	#return ac.accel_logarithmic(x,y,z)
	#return accel_halo_dm(x,y,z)
	#return accel_disk(x,y,z)


