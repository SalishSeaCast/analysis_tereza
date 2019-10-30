##calculates surface dic, no3, salinity for checking spinup

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


#changeable - at BR right now
ncname = 'PIACBC_2015_to_1127_SPINUP_massbal_comparison.nc'
start1 = '2015-01-01'
end1 = '2015-11-26'
st = dt.datetime(2015,1,1)
en = dt.datetime(2015,11,26)

sdir_1 = '/data/tjarniko/results/BASERUN_EXP/MAIN/PI_ACBC_2015/ncs/'
sdir_2 = '/data/tjarniko/results/BASERUN_EXP/MAIN/PI_ACBC_2015_2/ncs/'

y_st = st.timetuple().tm_yday
print(y_st)
y_en = en.timetuple().tm_yday
print(y_en)
ts_1st = np.arange(y_st,y_en+1,1)
days_in = len(ts_1st)

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
        nc_sens = sdir + '/SKOG_1d_*'+ ftype +'_T_' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        
    return sens_ar

sens_ar_1st = make_nclen(start1,end1,'carp', sdir_1)
sens_ar_2nd = make_nclen(start1,end1,'carp', sdir_2)
sens_ar_1st_grid = make_nclen(start1,end1,'grid', sdir_1)
sens_ar_2nd_grid = make_nclen(start1,end1,'grid', sdir_2)
sens_ar_1st_ptrc = make_nclen(start1,end1,'ptrc', sdir_1)
sens_ar_2nd_ptrc = make_nclen(start1,end1,'ptrc', sdir_2)

grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
vdir = grid['e2t'][0,:,:]
udir = grid['e1t'][0,:,:]
t_surfa = udir*vdir

def combine_files_surfa(files, surfa, var):
    stor_mol = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables[var][0,0,:,:]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * surfa
        totdic = np.sum(np.sum(var_tmp2))
        totdic_mols = totdic * (1/1000)
#         print(totflux)
        
        stor_mol[i] = totdic_mols
        i = i+1
    return stor_mol

surf_dic_1st = combine_files_surfa(sens_ar_1st,t_surfa, 'dissolved_inorganic_carbon')
surf_dic_2nd = combine_files_surfa(sens_ar_2nd,t_surfa, 'dissolved_inorganic_carbon')
surf_sal_1st = combine_files_surfa(sens_ar_1st_grid,t_surfa, 'vosaline')
surf_sal_2nd = combine_files_surfa(sens_ar_2nd_grid,t_surfa, 'vosaline')
surf_nit_1st = combine_files_surfa(sens_ar_1st_ptrc,t_surfa, 'nitrate')
surf_nit_2nd = combine_files_surfa(sens_ar_2nd_ptrc,t_surfa, 'nitrate')


f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
g.createDimension('days', days_in)
ts = g.createVariable('surf_dic_1st','f4',('days'))
ts[:] = surf_dic_1st
ts2 = g.createVariable('surf_dic_2nd','f4',('days'))
ts2[:] = surf_dic_2nd
ts3 = g.createVariable('surf_sal_1st','f4',('days'))
ts3[:] = surf_sal_1st
ts4 = g.createVariable('surf_sal_2nd','f4',('days'))
ts4[:] = surf_sal_2nd
ts5 = g.createVariable('surf_nit_1st','f4',('days'))
ts5[:] = surf_nit_1st
ts6 = g.createVariable('surf_nit_2nd','f4',('days'))
ts6[:] = surf_nit_2nd
f.close()
