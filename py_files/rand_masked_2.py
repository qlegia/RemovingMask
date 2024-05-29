# reading the Gaussian random field and mask it with kapp3_mask_Nside2048.fits
import numpy as np
import healpy as hp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import math
import scipy.io as sio

from colormap import cmbcmap
pi = math.pi

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

import sys

# set the parameters
# usage:  python rand_masked_2.py inst_num
inst_num = int(sys.argv[1]) # instance number from command line 
map_type = 'Linear'
prv_txt = 'masked'
Nside = 2048

dir = '../linearRF/'
fname = dir+map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + '.fits'
#rand_map = hp.read_map("randfield_Nside2048_instance1.fits")
rand_map = hp.read_map(fname)
mask = hp.read_map("kapp3_Mask_Nside2048.fits")

#mask = hp.read_map("kapp3_Mask_Nside2048.fits").astype(np.bool_)
#rand_map_masked = hp.ma(rand_map)
#rand_map_masked.mask = np.logical_not(mask)

rand_map_masked = np.multiply(rand_map, mask)

Lmax = 1200  
#computed the Fourier coefficients of the masked random field
alm = hp.map2alm(rand_map_masked,lmax=Lmax)

# save coefficients
sv = map_type + '_' + prv_txt + '_Nside' + str(Nside) + '_instance' + str(inst_num) + '.mat'
sio.savemat(sv, mdict={'alm':alm})

#plot figure
plt.figure(1)
cm = cmbcmap()

#sv_fig = "randfield_Nside2048_mask.png"
sv_fig = map_type + '_Nside' + str(Nside) + '_masked_instance'+str(inst_num) + '.png'
ti = "Masked Gaussian random field, instance " + str(inst_num)
#hp.mollview(rand_map_masked.filled(),title = ti, cmap=cm, min=-300, max=300, xsize=1200, nest=False)
hp.mollview(rand_map_masked,title = ti, cmap=cm, min=-100, max=100, xsize=1200, nest=False)
plt.title(ti)
plt.savefig(sv_fig,format='png',dpi=600)

