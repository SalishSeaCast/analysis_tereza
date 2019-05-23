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


#year BR
start = '2015-01-01'
end = '2015-12-31'

st = dt.datetime(2015,1,1)
en = dt.datetime(2015,12,31)

y_st = st.timetuple().tm_yday
print(y_st)
y_en = en.timetuple().tm_yday
print(y_en)
ts_BR = np.arange(y_st,y_en+1,1)

#PI year
start3 = '2015-01-01'
end3 = '2015-12-31'

st3 = dt.datetime(2015,1,1)
en3 = dt.datetime(2015,12,31)
y_st3 = st3.timetuple().tm_yday
print(y_st3)
y_en3 = en3.timetuple().tm_yday
print(y_en3)
ts_PI = np.arange(y_st3,y_en3+1,1)


sdir = '/data/tjarniko/results/BR_1st_2015/ncs/'
sdir3 = '/data/tjarniko/results/PREIND_1st_2015/ncs/'

def make_nclen(start,end,ftype, sdir):
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
        nc_sens = sdir + '/SKOG_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        sens_ar.append(tnc_sens[0])
    return sens_ar


BR_ar = make_nclen(start,end,'ptrc', sdir)
PI_ar = make_nclen(start3,end3,'ptrc', sdir3)


grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
vdir = grid['e2t'][0,0:878,20:398]
udir = grid['e1t'][0,0:878,20:398]
wdir = grid['e3t_0'][0,:,0:878,20:398]
wdir_20 = grid['e3t_0'][0,1:20,0:878,20:398]
surfa = vdir*udir
size_domain = wdir *surfa
size_domain_20 = wdir_20 *surfa

def calculate_total_N(files, size_domain):
    stor_nit = np.zeros(len(files))
    stor_am = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['nitrate'][0,:,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * size_domain
        totnit = np.sum(np.sum(var_tmp2))
        totnit_mols = totnit * (1/1000)        
        stor_nit[i] = totnit_mols
        
        var_tmp = G.variables['ammonium'][0,:,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * size_domain
        totam = np.sum(np.sum(var_tmp2))
        totam_mols = totam * (1/1000)        
        stor_am[i] = totam_mols
        i = i+1
    return stor_nit, stor_am


print('BR')
stor_nit_BR, stor_am_BR = calculate_total_N(BR_ar, size_domain)
stor_nit_PI, stor_am_PI = calculate_total_N(PI_ar, size_domain)

ncname = 'BRandPI2015_1styr_N_massbal.nc'
f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(NO3_mod))
g.createDimension('days', 365)
ts = g.createVariable('stor_nit_BR','f4',('days'))
ts[:] = stor_nit_BR
ts2 = g.createVariable('stor_am_BR','f4',('days'))
ts2[:] = stor_am_BR
ts3 = g.createVariable('stor_nit_PI','f4',('days'))
ts3[:] = stor_nit_PI
ts4 = g.createVariable('stor_am_PI','f4',('days'))
ts4[:] = stor_am_PI
f.close()



print('done')


