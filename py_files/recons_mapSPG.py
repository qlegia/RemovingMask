# load the coefficients from solving the regularization problem
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

Nside = 2048
Lmax = 100
# select approach 
fname = 'L100_SPG_noise.mat'
re_const_dat = sio.loadmat(fname)
# information about noise
#sig0 = re_const_dat['sigma0']
#kap = re_const_dat['kap']
# get the original values alm
org_alm = re_const_dat['org_alm']
print(org_alm.shape)
len_o = org_alm.shape[0]
org_alm = np.reshape(org_alm,(len_o,))
org_field = hp.alm2map(alms=org_alm, nside=Nside)

# get the reconstructed values alm
rec_alm = re_const_dat['rec_alm']
print(rec_alm.shape)
len_r = rec_alm.shape[0]
rec_alm = np.reshape(rec_alm,(len_r,))
reconst_field = hp.alm2map(alms=rec_alm, nside=Nside)

# get the power of fac = 10^{-pow}
re_pow1 = re_const_dat['pow']
re_pow = re_pow1[0][0]

# compute the error field
err_field = reconst_field - org_field

# compute the error
Npix = hp.nside2npix(Nside)
print(' Npix = ', Npix)
print('power=',re_pow, ' RMS error(hata) = ',np.linalg.norm(err_field)*2.0*np.sqrt(pi)/np.sqrt(Npix))
print('power=',re_pow,' RMS error(0) = ',np.linalg.norm(org_field)*2.0*np.sqrt(pi)/np.sqrt(Npix))
print('power=',re_pow,' relative RMS error = ',np.linalg.norm(err_field)/np.linalg.norm(org_field))

# compute the local errors
#load the mask
mask = hp.read_map("kapp3_Mask_Nside2048.fits")
idx0 = np.where(mask==0)  
idx1 = np.where(mask>0)  # same as np.where(a)
err0 = err_field[idx0]
err1 = err_field[idx1]
org0 = org_field[idx0]
org1 = org_field[idx1]
#
err0_hata = np.linalg.norm(err0)*2.0*np.sqrt(pi)/np.sqrt(len(err0))
err0_a0 = np.linalg.norm(org0)*2.0*np.sqrt(pi)/np.sqrt(len(org0))
rel0 = err0_hata/err0_a0
print('Number of points where mask=0', len(err0))
print("RMS Error where mask=0, err0(hata) = %2.3f" % err0_hata)
print('RMS Error where mask=0, err0(0) = %2.3f ' % err0_a0)
print('Relative err where mask=0 =      %2.3f  ' % rel0)
#
err1_hata = np.linalg.norm(err1)*2.0*np.sqrt(pi)/np.sqrt(len(err1))
err1_a0 = np.linalg.norm(org1)*2.0*np.sqrt(pi)/np.sqrt(len(org1))
rel1 = err1_hata/ err1_a0
print('Number of points where mask>0', len(err1))
print('RMS Error where mask=1, err0(hata) = %2.3f' % err1_hata)
print('RMS Error where mask=1, err0(0) = %2.3f' % err1_a0)
print('Relative err where mas>=0 =  %2.3f'  % rel1)

# saving and plotting
#ti = r'Reconstructed field with no noise'
#sv_fig = 'reconstructed_Gaussian_field_no_noise_QR.png'
ti = r'Reconstructed field with Gaussian noise with power spectrum $= 10^{-'+str(re_pow)+'} C_{\ell}$'
sv_fig = 'reconstr_field_SPG_gauss_noise_1e-'+str(re_pow) + '.png'

cm = cmbcmap()
hp.mollview(reconst_field, title = ti, cmap=cm, min=-70, max=70, xsize=1200, nest=False)
plt.title(ti)
plt.savefig(sv_fig,format='png',dpi=600) 

#err_ti = 'Point-wise error field with no noise'
#err_sv = 'error_field_no_noise_QR.png'
err_ti = 'Point-wise error field'
err_sv = 'error_field_SPG_gauss_noise_1e-'+str(re_pow) +'.png'

cm = cmbcmap()
hp.mollview(err_field, title = err_ti, cmap=cm, min=-70, max=70, xsize=1200, nest=False)
plt.title(err_ti)
plt.savefig(err_sv,format='png',dpi=600) 

#sv_fits = 'reconstructed_field'+str(appr) + '.fits'
# save into a fits file as well
#hp.write_map(sv_fits,re_const_field,overwrite=True,dtype='float64')


