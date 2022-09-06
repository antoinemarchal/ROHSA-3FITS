import os 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from ROHSApy import ROHSA, fit_spec

import marchalib as ml

#Open data
#Open data
path = "../data/"
fitsname = "GHIGLS_DFN_Tb_32.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
cube = hdu[0].data

fitsname = "GHIGLS_DFN_Tb_32_rms.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
rms = hdu[0].data

n_gauss = 3
lb_sig_init = 1.
ub_sig_init = 10.                                                                                                       
iprint_init = 1                                                 
amp_fact_init = 0.66                                                                                                      
sig_init = 2.                                                                                                       
maxiter_init = 15000                                                                                              

vmin = [0]    
vmax = [cube.shape[0]]    

Tb = np.nanmean(cube,(1,2))
core_init = fit_spec(Tb=Tb)
rms_init = 1.# np.sqrt(np.sum(rms**2, (1,2))) / (rms.shape[2]*rms.shape[1])
init_array = core_init.init_spectrum(np.full((3*n_gauss),1.), n_gauss, Tb, lb_sig_init, ub_sig_init, 
                            iprint_init, amp_fact_init, sig_init, maxiter_init, rms_init, vmin, vmax)    
model = core_init.model_spectrum(init_array, Tb, n_gauss)



