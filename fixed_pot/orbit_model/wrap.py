#Wrappers
import numpy as np
import config as cfg


def phi_2_wrap(x):
    array = np.zeros(len(x))
    for i in range(0,len(x)):
        if(x[i]>cfg.phi_1_min and x[i]<cfg.phi_1_max ):
            array[i]= cfg.phi_2_spl(x[i])
        else:
            array[i]= np.Inf
    return array
        
def d_hel_wrap(x):
    array = np.zeros(len(x))
    for i in range(0,len(x)):
        if(x[i]>cfg.phi_1_min and x[i]<cfg.phi_1_max ):
            array[i]= cfg.d_hel_spl(x[i])
        else:
            array[i]= np.Inf
    return array
    
def v_hel_wrap(x):
    array = np.zeros(len(x))
    for i in range(0,len(x)):
        if(x[i]>cfg.phi_1_min and x[i]<cfg.phi_1_max ):
            array[i]= cfg.v_hel_spl(x[i])
        else:
            array[i]= np.Inf
    return array
        
def mu_1_wrap(x):
    array = np.zeros(len(x))
    for i in range(0,len(x)):
        if(x[i]>cfg.phi_1_min and x[i]<cfg.phi_1_max ):
            array[i]= cfg.mu_1_spl(x[i])
        else:
            array[i]= np.Inf
    return array

def mu_2_wrap(x):
    array = np.zeros(len(x))
    for i in range(0,len(x)):
        if(x[i]>cfg.phi_1_min and x[i]<cfg.phi_1_max ):
            array[i]= cfg.mu_2_spl(x[i])
        else:
            array[i]= np.Inf
    return array

def mu_ra_wrap(x):
    array = np.zeros(len(x))
    for i in range(0,len(x)):
        if(x[i]>cfg.phi_1_min and x[i]<cfg.phi_1_max ):
            array[i]= cfg.mu_ra_spl(x[i])
        else:
            array[i]= np.Inf
    return array

def mu_dec_wrap(x):
    array = np.zeros(len(x))
    for i in range(0,len(x)):
        if(x[i]>cfg.phi_1_min and x[i]<cfg.phi_1_max ):
            array[i]= cfg.mu_dec_spl(x[i])
        else:
            array[i]= np.Inf
    return array
