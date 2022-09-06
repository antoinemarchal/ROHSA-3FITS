import os 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from ROHSApy import ROHSA, fit_spec

import marchalib as ml

#Open data
#Open data
path = "../data/"
fitsname = "GHIGLS_DFN_Tb.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
cube = hdu[0].data[0][100:400,:32,:32]
rms2d = np.std(cube[:20], 0)

rms = np.array([rms2d * (1+channel/20) for channel in cube])

filename_parameters = path + "PARAMETERS/GHIGLS_DFN_Tb_parameters_run_0_ifort.txt"
filename = path + "DAT/GHIGLS_DFN_Tb.dat"
fileout = path + "ROHSA/GHIGLS_DFN_Tb_gauss_run_0_ifort.dat"
filename_noise = path + "DAT/GHIGLS_DFN_Tb_rms.dat"
n_gauss = 3
lambda_amp = 1.                                                                                                                   
lambda_mu = 1.
lambda_sig = 1.                                                                              
lambda_var_sig = 1. 
amp_fact_init = 0.66                                                                                                      
sig_init = 2.                                                                                                       
lb_sig_init = 1.
ub_sig_init = 10.                                                                                                       
lb_sig = 1.
ub_sig = 100.                                                                                                       
init_option = 'mean'                                                                                            
maxiter_init = 15000                                                                                              
maxiter = 800
m = 10                                                                                                
noise = ".true."                                                                                                              
regul = ".true."                                                                      
descent = ".true."                                                                                                                 
lstd = 1                                                                                                                        
ustd = 1                                                                                                                       
iprint = -1                                                                                                                     
iprint_init = -1                                                 
save_grid = ".false."  
init_spec = ".false."

core = ROHSA(cube, hdr=hdr)
core_rms = ROHSA(rms)

core.cube2dat(filename=filename)
core_rms.cube2dat(filename=filename_noise)
core.gen_parameters_3D(filename_parameters=filename_parameters,
                    filename=filename, 
                    fileout=fileout,  
                    filename_noise = filename_noise,
                    # filename_init_spec = filename_init_spec,
                    n_gauss=n_gauss,
                    amp_fact_init=amp_fact_init,
                    lambda_amp=lambda_amp,
                    lambda_mu=lambda_mu,
                    lambda_sig=lambda_sig,
                    lambda_var_sig=lambda_var_sig,
                    maxiter=maxiter,
                    sig_init=sig_init,
                    lb_sig=lb_sig,
                    ub_sig=ub_sig,
                    lb_sig_init=lb_sig_init,
                    ub_sig_init=ub_sig_init,
                    init_option=init_option,
                    noise=noise,
                    descent=descent,
                    lstd=lstd,
                    ustd=ustd,
                    iprint=iprint,
                    iprint_init=iprint_init,
                    save_grid=save_grid,
                    init_spec=init_spec)
