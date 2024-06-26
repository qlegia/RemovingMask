# Usage coeff_rand_field_v5.py instance
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import math
import scipy.io as sio
#import os
import sys
from colormap import cmbcmap
pi = math.pi

subdir = '/'
map_type = 'Noise_1e_1'
map_txt = 'Gaussian noise'
Nside = 2048
L = 100
halfL = 50
prv_txt = 'Lmax100'

fac = 1e-1
cl_bf = fac*np.ones([L+1])

for el in range(halfL+1,L+1):
   #cl_bf[el] = -2*el/L + 2
   cl_bf[el] = fac*(-2*el/(L+1) + 2)

T_ran = hp.synfast(cl_bf,Nside)
# Fourier coefficients of random field
alm = hp.map2alm(T_ran,lmax=L)

#%% save coefficients
inst_num = int(sys.argv[1]) # instance number from command line 
sv = map_type + '_' + 'Nside' + str(Nside) + '_instance' + str(inst_num) +'.mat'
sio.savemat(sv, mdict={'alm':alm})

# random field from alm
randfield = hp.alm2map(alms=alm,nside=Nside)
print('min/max randfield', min(randfield), max(randfield))
#%% plot figure
# random field from alm
sv_fig = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + '.png'
sv_fits = map_type + '_Nside' + str(Nside) + '_instance' + str(inst_num) + '.fits'
plt.figure(1)
cm = cmbcmap()
ti = map_txt + r' with angular power spectrum $\Upsilon_{\ell} = 10^{-2} C_{\ell}$, Nside $= %d$' % Nside
#ti = map_txt + ' with angular power spectrum $\Upsilon_{\ell} = 10^{-2} C_{\ell}$, Nside $= 2048$'
hp.mollview(randfield, title = ti, cmap=cm, min=-6, max=6, xsize=1200, nest=False)
plt.title(ti)
plt.savefig(sv_fig,format='png',dpi=600)

# save into a fits file as well
hp.write_map(sv_fits,randfield,overwrite=True,dtype='float64')



