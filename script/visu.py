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

core = ROHSA(cube, hdr=hdr)

filename = "ROHSA/GHIGLS_DFN_Tb_gauss_run_0_gfort.dat"
output = core.read_gaussian(path+filename)    

fitsname = "ROHSA/GHIGLS_DFN_Tb_gauss_run_0_ifort.fits"
hdu = fits.open(path+fitsname)
hdr = hdu[0].header
output2 = hdu[0].data



