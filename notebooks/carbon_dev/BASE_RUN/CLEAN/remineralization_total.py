from __future__ import print_function
from numpy import *
from scipy import *
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use('seaborn-whitegrid')
import netCDF4 as nc

import cmocean as cm
import glob
import sys
sys.path.append('/data/tjarniko/mocsy')
sys.path.append('/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/CCCmaDEV/CCCma_src')
import mocsy
import CCCma
import CCCma_stations as cs
from matplotlib import reload
import arrow
import gsw


def make_nclen(start,end,ftype):
    base_ar = []
    sens_ar = []
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    dayslen = len(arrow_array)
    for i in range(0,dayslen):
        tdate = arrow_array[i][0]
        ddmmmyy = tdate.format('DDMMMYY').lower()
        ymd = tdate.format('YYYYMMDD')
        nc_sens = '/results/SalishSea/hindcast.201812/' + ddmmmyy + '/SalishSea_1h_*'+ftype+'*.nc'
        tnc_sens = glob.glob(nc_sens)
        #print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        
    return sens_ar

start = '2015-01-01'
end = '2015-12-31'
grid_ar = make_nclen(start,end,'grid_T')
prod_ar = make_nclen(start,end,'prod_T')
ptrc_ar = make_nclen(start,end,'ptrc_T')

grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
vdir = grid['e2t'][0,0:878,20:398]
udir = grid['e1t'][0,0:878,20:398]
wdir = grid['e3t_0'][0,:,0:878,20:398]
surfa = vdir*udir
size_domain = wdir *surfa


import timeit
start_time = timeit.default_timer()
print(start_time)

year_remin = np.zeros(365)
rate_remin = 0.1987
for i in range(0,365):
    print(i)
    ptrc_test = nc.Dataset(ptrc_ar[i])
    grid_test = nc.Dataset(grid_ar[i])
    DON = ptrc_test['dissolved_organic_nitrogen'][:]
    PON = ptrc_test['particulate_organic_nitrogen'][:]
    temp = grid_test['votemper'][:]
    temp_dayavg = np.mean(temp, axis = 0)
    tgfunc = np.exp(0.07 * (temp_dayavg-20) )
    DON_dayavg = np.mean(DON, axis = 0)
    PON_dayavg = np.mean(PON, axis = 0)
    DON_daily_remin = DON_dayavg[:,0:878,20:398] * size_domain * rate_remin * 1/1000 * 6.625 * tgfunc[:,0:878,20:398]
    PON_daily_remin = PON_dayavg[:,0:878,20:398] * size_domain * rate_remin * 1/1000 * 6.625 * tgfunc[:,0:878,20:398]

    total_daily_remin = (np.sum(np.sum(np.sum(DON_daily_remin)))) + (np.sum(np.sum(np.sum(PON_daily_remin))))
    year_remin[i] = total_daily_remin
    
end_time = timeit.default_timer()
print(end_time - start_time)

f = nc.Dataset('REM_2015.nc','w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(NO3_mod))
g.createDimension('days', 365)
ts = g.createVariable('year_remin','f4',('days'))
ts[:] = year_remin

f.close()