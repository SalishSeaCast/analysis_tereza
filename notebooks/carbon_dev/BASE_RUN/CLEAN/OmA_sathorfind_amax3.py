from __future__ import print_function
from numpy import *
from scipy import *
import netCDF4 as nc
import numpy as np
import scipy as sp
import seawater
import datetime as dt

""
from salishsea_tools import (
    nc_tools,
    viz_tools, 
    geo_tools,
    tidetools
)
from salishsea_tools import visualisations as vis
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import matplotlib.patches as patches
plt.style.use('seaborn-whitegrid')
import netCDF4 as nc

import cmocean as cm
import glob
import sys
sys.path.append('/data/tjarniko/mocsy')
sys.path.append('/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/CCCmaDEV/CCCma_src')
sys.path.append('/data/tjarniko/MEOPAR/tools/SalishSeaTools/salishsea_tools/')
import mocsy
import river_201702 as rv
import mocsy
import CCCma
import CCCma_stations as cs
from matplotlib import reload
import arrow
import gsw
import datetime as dt

import timeit
start_time = timeit.default_timer()
print(start_time)

OmA = nc.Dataset('OmA_2015.nc')
#switch BR and PI
BR_omA = OmA['model_output']['OmAr_pi']
PI_omA = OmA['model_output']['OmAr_br']

t_nc = nc.Dataset('/results2/SalishSea/nowcast-green.201806/01jan18/SalishSea_1h_20180101_20180101_grid_T.nc')
zlevels = (t_nc['deptht'][:])


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def find_oma_depths_via_amax_2(Omega_A_ar, zlevels):
    #subtract 1 to find values close to sat. horizon.
#     oma_close_to_1 = np.abs(Omega_A_ar[:,:,:] - 1)
    oma_depths = np.zeros([898,398])
    depth_inds = np.argmax(Omega_A_ar<1, axis=0)
    for i in range(0,898):
        for j in range(0,398):
            test_slice = Omega_A_ar[:,i,j]
            if np.min(test_slice) > 1:
                oma_depths[i,j] = np.nan
            else:
                oma_depths[i,j] = zlevels[depth_inds[i,j]]
            
    return oma_depths

def find_oma_depths_via_susan(Omega_A_ar,tmask):
    depthindex = np.argmax(Omega_A_ar<1, axis=0)
    maxsal = np.amax(Omega_A_ar, axis=0)<1
    t_depths = np.ma.array(depthindex, mask=(1 - (tmask[0] * maxsal)));
    oma_depths = np.zeros([898,398])
    for i in range(0,898):
        for j in range(0,398):
            if np.ma.is_masked(t_depths[i,j]):
                oma_depths[i,j] = np.nan
            else:
                oma_depths[i,j] = zlevels[t_depths[i,j]]
    
    return oma_depths

OmA = nc.Dataset('OmA_2015.nc')
BR_omA = OmA['model_output']['OmAr_pi']
PI_omA = OmA['model_output']['OmAr_br']

tmask = nc.Dataset('/home/sallen/MEOPAR/grid/mesh_mask201702.nc')['tmask'][0]
t_nc = nc.Dataset('/results2/SalishSea/nowcast-green.201806/01jan18/SalishSea_1h_20180101_20180101_grid_T.nc')
zlevels = (t_nc['deptht'][:])

oma_d_PI = np.zeros([365,898,398])
oma_d_BR = np.zeros([365,898,398])
for day in range(0,365):
    print(day)
    OmA_BR_test = BR_omA[day,:,:,:]
    OmA_PI_test = PI_omA[day,:,:,:]
    
    oma_dep_PI = find_oma_depths_via_susan(OmA_PI_test,tmask)
    oma_dep_BR = find_oma_depths_via_susan(OmA_BR_test,tmask)
    oma_d_PI[day,:,:] = oma_dep_PI
    oma_d_BR[day,:,:] = oma_dep_BR
    

f = nc.Dataset('OmA_horizon_2015_susanattempt.nc','w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(3))
g.createDimension('days', 365)
g.createDimension('ys', 898)
g.createDimension('xs', 398)
ts = g.createVariable('OmArHORIZON_pi','f4',('days','ys','xs'))
ts[:] = oma_d_PI
ts2 = g.createVariable('OmArHORIZON_br','f4',('days','ys','xs'))
ts2[:] = oma_d_BR
f.close()

elapsed = timeit.default_timer() - start_time
print(elapsed)