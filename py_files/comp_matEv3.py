# Usage: python comp_matEv3 L1max L2max m
# generate a matrix E of size (L1max+L2max+1-m, L1max+1-m)
import numpy as np
import healpy as hp
import time
from sympy.physics.wigner import gaunt
import sys

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

# to save matrix E as matlab .mat file
import scipy.io as sio

# read values from command line
L1max = int(sys.argv[1])
L2max = int(sys.argv[2]) 
m     = int(sys.argv[3])

# load vlm from pre-computed mat file
vfname = '../mat_files/v_ell0_500_gaussian.mat'
vLmax = 500
mat = sio.loadmat(vfname)
vl = np.reshape(mat['vl'],[vLmax+1])

# compute the matrix E
matE = np.zeros([L1max+L2max+1-m,L1max+1-m],dtype=float)
start = time.time()
# 
for ell in range(m,L1max+L2max+1):
    for ell_1 in range(m,L1max+1):
        if ( (ell+ell_1)%2 == 0):
            res =  0
            lower = np.abs(ell-ell_1)
            upper = np.min([L2max+1,ell+ell_1+1])
            #for ell_2 in range(0,L2max+1):
            for ell_2 in range(lower,upper):
               if ( ell_2 % 2 == 0):
               # so (ell + ell_1 + ell_2 even number
                   res = res + (-1)**m*gaunt(ell,ell_1,ell_2,-m,m,0)*vl[ell_2]
            matE[ell-m,ell_1-m] = res

# save matrix E as a .mat file
fname = 'E_L1max'+str(L1max) + '_L2max'+str(L2max)+'_m'+str(m)+'.mat'
sio.savemat(fname, mdict={'E':matE}) 
finish = time.time()
tt = finish - start
print("Elapsed time = %4.2f  "% tt)

