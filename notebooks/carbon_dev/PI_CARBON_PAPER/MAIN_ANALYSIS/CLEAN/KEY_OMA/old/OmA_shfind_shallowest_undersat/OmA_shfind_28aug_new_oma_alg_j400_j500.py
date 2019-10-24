from __future__ import print_function


j_st = 400
j_en = 500
ncname = 'OmA_horizon_2015_newalg_' + str(j_st) + '_' + str(j_en) + '.nc'
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

def find_depth(dp,prof):
    #finds saturation horizon given a profile and corresponding depths
    first_proper_undersat = np.nan
    depth_undersat = np.nan    
    dummy_var = 0
    #tot masked
    if prof.mask.all():
        dummy_var = 0
        depth_undersat = np.nan
        #print('masks all around!')
    elif np.ma.min(prof) >1:
        depth_undersat = np.nan
        #print('saturated column')
    elif np.ma.max(prof) <1:
        depth_undersat = 0
        #print('undersat to surface')        
    else:
        t_ind = np.where(prof<1)
        t_indar = t_ind[0][0]
        t_indss = np.where(prof>=1)
        t_indsssar = t_indss[0][0]
        if t_indar.size == 0:
            dummy_var = 0
        else:
            if (t_indar.size != 0) & (t_indsssar.size == 0):
                depth_undersat = 0
                first_proper_undersat = 0
                dummy_var = 0
                #print('undersat to surface!')
                max_supsat = np.nan
            else:    
                max_supsat = np.max(t_indsssar)    
                try:
                    first_proper_undersat = np.min(t_indar)
                except:
                    dummy_var = 0
                    #print("An exception occurred")
                if first_proper_undersat == 0:
                    depth_undersat = dp[0]
                if np.isnan(first_proper_undersat):
                    dummy_var = 0
                    #print('saturated watercolumn!')
                else:
                    depth_undersat = (dp[first_proper_undersat]+dp[first_proper_undersat-1])/2
    return depth_undersat

OmA = nc.Dataset('OmA_2015.nc')
#switch BR and PI
BR_omA = OmA['model_output']['OmAr_pi']
PI_omA = OmA['model_output']['OmAr_br']

t_nc = nc.Dataset('/results2/SalishSea/nowcast-green.201806/01jan18/SalishSea_1h_20180101_20180101_grid_T.nc')
zlevels = (t_nc['deptht'][:])


import timeit
start_time = timeit.default_timer()
print(start_time)

oma_d_PI = np.zeros([365,898,398])
oma_d_BR = np.zeros([365,898,398])
for day in range(0,365):
    print(day)
    for j in range(j_st,j_en):
        if j%10 == 0:
            print(j)
        for i in range(0,398):
            OmA_BR_test = BR_omA[day,:,j,i]
            OmA_PI_test = PI_omA[day,:,j,i]
            
            oma_dep_BR = find_depth(zlevels,OmA_BR_test)
            oma_dep_PI = find_depth(zlevels,OmA_PI_test)
            oma_d_PI[day,j,i] = oma_dep_PI
            oma_d_BR[day,j,i] = oma_dep_BR
            if (i == 250) & (j%10 == 0):
                print('BR depth: '+str(oma_dep_BR), 'PI depth: '+str(oma_dep_PI))       

elapsed = timeit.default_timer() - start_time
print(elapsed)
    

f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
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