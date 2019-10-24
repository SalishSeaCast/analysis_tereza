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

start = '2015-01-01'
end = '2015-12-31'

st = dt.datetime(2015,1,1)
en = dt.datetime(2015,12,31)


y_st = st.timetuple().tm_yday
print(y_st)
y_en = en.timetuple().tm_yday
print(y_en)


bdir = '/data/tjarniko/results/PREIND_1st_2015/ncs/'
sdir = '/data/tjarniko/results/BR_1st_2015/ncs'

thalweg_file='/home/sallen/MEOPAR/Tools/bathymetry/thalweg_working.txt'
thalweg_pts = np.loadtxt(thalweg_file, delimiter=' ', dtype=int)

figstring = 'BR2015'

def make_nclen(start,end,ftype, bdir, sdir):
    base_ar = []
    sens_ar = []
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)
    
    #print('days: '+str(len(arrow_array)))
    dayslen = len(arrow_array)
    hovmoller_base = np.zeros([1533,40,dayslen])
    hovmoller_sens = np.zeros([1533,40,dayslen])

    for i in range(0,dayslen):
        tdate = arrow_array[i][0]
        ddmmmyy = tdate.format('DDMMMYY').lower()
        ymd = tdate.format('YYYYMMDD')
        nc_base = bdir + '/SKOG_1d_*'+ ftype +'_T_' + ymd + '-' + ymd + '.nc'
        nc_sens = sdir + '/SKOG_1d_*'+ ftype +'_T_' + ymd + '-' + ymd + '.nc'
        tnc_base = glob.glob(nc_base) 
        tnc_sens = glob.glob(nc_sens)
        #print(tnc_sens)
        
        base_ar.append(tnc_base[0])
        sens_ar.append(tnc_sens[0])
        
    return base_ar, sens_ar

base_ar, sens_ar = make_nclen(start,end,'carp', bdir, sdir)
base_arG, sens_arG = make_nclen(start,end,'grid', bdir, sdir)

def combine_files(files, var, jss, iss):

    time = np.array([])
    var_list = []
    i = 0
    for f in files:
        if i%30 == 0:
            print(i)
        G = nc.Dataset(f)
#         if i == 1:
#             print(G)
        var_tmp = G.variables[var][:]
        if var == 'co2_flux_mmol_m2_s':
            var_tmp=var_tmp[:,jss,iss]
        
        else:
            var_tmp=var_tmp[:,:,jss,iss]
        
        var_list.append(var_tmp)
        i = i+1

    var_ary = np.concatenate(var_list, axis=0)
    return var_ary, time

dic_base = combine_files(base_ar,'dissolved_inorganic_carbon',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
dic_base2 = dic_base[0]
print(np.shape(dic_base2))

flux_base = combine_files(base_ar,'co2_flux_mmol_m2_s',thalweg_pts[:, 0],thalweg_pts[:, 1])
flux_base2 = flux_base[0]
print(np.shape(flux_base2))

temp_base = combine_files(base_arG,'votemper',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
temp_base2 = temp_base[0]
print(np.shape(temp_base2))

sal_base = combine_files(base_arG,'vosaline',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
sal_base2 = sal_base[0]
print(np.shape(sal_base2))

ta_base = combine_files(base_ar,'total_alkalinity',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
ta_base2 = ta_base[0]

dic_sens = combine_files(sens_ar,'dissolved_inorganic_carbon',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
dic_sens2 = dic_sens[0]
print(np.shape(dic_sens2))

flux_sens = combine_files(sens_ar,'co2_flux_mmol_m2_s',thalweg_pts[:, 0],thalweg_pts[:, 1])
flux_sens2 = flux_sens[0]
print(np.shape(flux_sens2))

temp_sens = combine_files(sens_arG,'votemper',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
temp_sens2 = temp_sens[0]
print(np.shape(temp_sens2))

sal_sens = combine_files(sens_arG,'vosaline',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
sal_sens2 = sal_sens[0]
print(np.shape(sal_sens2))

ta_sens = combine_files(sens_ar,'total_alkalinity',thalweg_pts[:, 0],thalweg_pts[:, 1])
#extract masked array from tuple
ta_sens2 = ta_sens[0]

ncname = 'BR2015_1styr_hovmoller.nc'
f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(NO3_mod))
g.createDimension('days', 365)
g.createDimension('depths', 40)
g.createDimension('thwg_pt', 1533)
ts = g.createVariable('dic_hovmol','f4',('days','depths','thwg_pt'))
ts[:] = dic_sens2
ts2 = g.createVariable('ta_hovmol','f4',('days','depths','thwg_pt'))
ts2[:] = ta_sens2
ts2 = g.createVariable('temp_hovmol','f4',('days','depths','thwg_pt'))
ts2[:] = temp_sens2
ts2 = g.createVariable('sal_hovmol','f4',('days','depths','thwg_pt'))
ts2[:] = sal_sens2
ts2 = g.createVariable('flux_hovmol','f4',('days','thwg_pt'))
ts2[:] = flux_sens2
f.close()

ncname = 'PI2015_1styr_hovmoller.nc'
f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(NO3_mod))
g.createDimension('days', 365)
g.createDimension('depths', 40)
g.createDimension('thwg_pt', 1533)
ts = g.createVariable('dic_hovmol','f4',('days','depths','thwg_pt'))
ts[:] = dic_base2
ts2 = g.createVariable('ta_hovmol','f4',('days','depths','thwg_pt'))
ts2[:] = ta_base2
ts2 = g.createVariable('temp_hovmol','f4',('days','depths','thwg_pt'))
ts2[:] = temp_base2
ts2 = g.createVariable('sal_hovmol','f4',('days','depths','thwg_pt'))
ts2[:] = sal_base2
ts2 = g.createVariable('flux_hovmol','f4',('days','thwg_pt'))
ts2[:] = flux_base2
f.close()