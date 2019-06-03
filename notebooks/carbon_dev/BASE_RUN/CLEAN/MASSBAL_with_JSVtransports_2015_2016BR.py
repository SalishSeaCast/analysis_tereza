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
end = '2015-12-30'
start3 = '2016-01-01'
end3 = '2016-12-30'

st = dt.datetime(2015,1,1)
en = dt.datetime(2015,12,30)
st3 = dt.datetime(2016,1,1)
en3 = dt.datetime(2016,12,30)

ncname_BR = 'MASSBAL_BR2015_to1230_spunup.nc'
ncname_BR2 = 'MASSBAL_BR2016_to1230_spunup.nc'

sdir = '/data/tjarniko/results/BR_2nd_2015_cop/SKOG_2/ncs/'
sdir2 = '/data/tjarniko/results/BR_2016/ncs/'

y_st = st.timetuple().tm_yday
print(y_st)
y_en = en.timetuple().tm_yday
print(y_en)
ts_BR = np.arange(y_st,y_en+1,1)
days_in_2015 = np.size(ts_BR)

#BR2 year
y_st3 = st3.timetuple().tm_yday
print(y_st3)
y_en3 = en3.timetuple().tm_yday
print(y_en3)
ts_BR2= np.arange(y_st3,y_en3+1,1)
days_in_2016 = np.size(ts_BR2)




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
        print(tnc_sens[0])
    return sens_ar

def make_nclen_transport(start,end,sdir):
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
        nc_sens = sdir + '/SKOG_1d_*' +'dian_U_' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        sens_ar.append(tnc_sens[0])

    return sens_ar

def make_nclen_transport_V(start,end,sdir):
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
        nc_sens = sdir + '/SKOG_1d_*' +'dian_V_' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        sens_ar.append(tnc_sens[0])

    return sens_ar

BR_ar = make_nclen(start,end,'carp', sdir)
BR_ar_tp = make_nclen_transport(start,end, sdir)
BR_ar_tpV = make_nclen_transport_V(start,end, sdir)

BR_ar2 = make_nclen(start3,end3,'carp', sdir2)
BR_ar2_tp = make_nclen_transport(start3,end3, sdir2)
BR_ar2_tpV = make_nclen_transport_V(start3,end3, sdir2)

grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
vdir = grid['e2t'][0,0:878,20:398]
udir = grid['e1t'][0,0:878,20:398]
wdir = grid['e3t_0'][0,:,0:878,20:398]
wdir_20 = grid['e3t_0'][0,0:20,0:878,20:398]
wdir_20_100 = grid['e3t_0'][0,20:27,0:878,20:398]
wdir_deep = grid['e3t_0'][0,27:40,0:878,20:398]
surfa = vdir*udir
size_domain = wdir *surfa
size_domain_20 = wdir_20 *surfa
size_domain_20_100 = wdir_20_100 *surfa
size_domain_deep = wdir_deep *surfa

def calculate_total_C(files, size_domain):
    stor_mol = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        print(f)
        var_tmp = G.variables['dissolved_inorganic_carbon'][0,:,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * size_domain
        totdic = np.sum(np.sum(var_tmp2))
        totdic_mols = totdic * (1/1000)        
        stor_mol[i] = totdic_mols
        i = i+1
    return stor_mol

def calculate_surface_C(files, surfa):
    stor_mol = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['dissolved_inorganic_carbon'][0,0,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * surfa
        totdic = np.sum(np.sum(var_tmp2))
        totdic_mols = totdic * (1/1000)        
        stor_mol[i] = totdic_mols
        i = i+1
    return stor_mol

def calculate_surface20_C(files, surfa):
    stor_mol = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['dissolved_inorganic_carbon'][0,0:20,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * size_domain_20
        totdic = np.sum(np.sum(var_tmp2))
        totdic_mols = totdic * (1/1000)        
        stor_mol[i] = totdic_mols
        i = i+1
    return stor_mol

def calculate_20_100_C(files, surfa):
    stor_mol = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['dissolved_inorganic_carbon'][0,20:27,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * size_domain_20_100
        totdic = np.sum(np.sum(var_tmp2))
        totdic_mols = totdic * (1/1000)        
        stor_mol[i] = totdic_mols
        i = i+1
    return stor_mol

def calculate_100_deep_C(files, surfa):
    stor_mol = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['dissolved_inorganic_carbon'][0,27:40,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * size_domain_deep
        totdic = np.sum(np.sum(var_tmp2))
        totdic_mols = totdic * (1/1000)        
        stor_mol[i] = totdic_mols
        i = i+1
    return stor_mol

def calculate_flux(files, surfa):
    stor_flx = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['co2_flux_mmol_m2_s'][0,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * surfa
        totflux = np.sum(np.sum(var_tmp2))
        #mol/day
        totflux_daily_moles = totflux * 60 * 60 * 24 * (1/1000)
        stor_flx[i] = totflux_daily_moles
        i = i+1

    return stor_flx

def calculate_transports(files):
    stor_trans = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['DIC_UT'][:,:,:,20]
        var_tmp[var_tmp == 1e+20] = 0
        #mmol/s > mol/day
        var_tmp2 = np.sum(var_tmp)*(1/1000)*60*60*24
        stor_trans[i] = var_tmp2
        i = i+1

    return stor_trans

def calculate_transports_JS(files):
    stor_trans_JS = np.zeros(len(files))

    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
        var_tmp = G.variables['DIC_VT'][:,:,878,0:120]
        var_tmp[var_tmp == 1e+20] = 0
        #mmol/s > mol/day
        var_tmp2 = np.sum(var_tmp)*(1/1000)*60*60*24
        stor_trans_JS[i] = var_tmp2
        i = i+1

    return stor_trans_JS

print('BR')
stor_mol_BR = calculate_total_C(BR_ar, size_domain)
stor_mol_surf_BR = calculate_surface_C(BR_ar, surfa)
stor_mol_20_BR = calculate_surface20_C(BR_ar, size_domain_20)
stor_mol_20_100_BR = calculate_20_100_C(BR_ar, size_domain_20_100)
stor_mol_deep_BR = calculate_100_deep_C(BR_ar, size_domain_deep)
stor_flx_BR = calculate_flux(BR_ar, surfa)
stor_trans_BR = calculate_transports(BR_ar_tp)
stor_trans_BR_JS = calculate_transports_JS(BR_ar_tpV)


f = nc.Dataset(ncname_BR,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(NO3_mod))
g.createDimension('days', days_in_2015)
ts = g.createVariable('stor_mol_BR','f4',('days'))
ts[:] = stor_mol_BR
ts2 = g.createVariable('stor_mol_surf_BR','f4',('days'))
ts2[:] = stor_mol_surf_BR
ts3 = g.createVariable('stor_mol_20_BR','f4',('days'))
ts3[:] = stor_mol_20_BR
ts3b = g.createVariable('stor_mol_20_100_BR','f4',('days'))
ts3b[:] = stor_mol_20_100_BR
ts3c = g.createVariable('stor_mol_deep_BR','f4',('days'))
ts3c[:] = stor_mol_deep_BR
ts4 = g.createVariable('stor_flx_BR','f4',('days'))
ts4[:] = stor_flx_BR
ts5 = g.createVariable('stor_trans_BR','f4',('days'))
ts5[:] = stor_trans_BR
ts6 = g.createVariable('stor_trans_BR_JS','f4',('days'))
ts6[:] = stor_trans_BR_JS
f.close()

print('BR2')
stor_mol_BR = calculate_total_C(BR_ar2, size_domain)
stor_mol_surf_BR = calculate_surface_C(BR_ar2, surfa)
stor_mol_20_BR = calculate_surface20_C(BR_ar2, size_domain_20)
stor_mol_20_100_BR = calculate_20_100_C(BR_ar2, size_domain_20_100)
stor_mol_deep_BR = calculate_100_deep_C(BR_ar2, size_domain_deep)
stor_flx_BR = calculate_flux(BR_ar2, surfa)
stor_trans_BR = calculate_transports(BR_ar2_tp)
stor_trans_BR_JS = calculate_transports_JS(BR_ar2_tpV)


f = nc.Dataset(ncname_BR2,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(NO3_mod))
g.createDimension('days', days_in_2016)
ts = g.createVariable('stor_mol_BR','f4',('days'))
ts[:] = stor_mol_BR
ts2 = g.createVariable('stor_mol_surf_BR','f4',('days'))
ts2[:] = stor_mol_surf_BR
ts3 = g.createVariable('stor_mol_20_BR','f4',('days'))
ts3[:] = stor_mol_20_BR
ts3b = g.createVariable('stor_mol_20_100_BR','f4',('days'))
ts3b[:] = stor_mol_20_100_BR
ts3c = g.createVariable('stor_mol_deep_BR','f4',('days'))
ts3c[:] = stor_mol_deep_BR
ts4 = g.createVariable('stor_flx_BR','f4',('days'))
ts4[:] = stor_flx_BR
ts5 = g.createVariable('stor_trans_BR','f4',('days'))
ts5[:] = stor_trans_BR
ts6 = g.createVariable('stor_trans_BR_JS','f4',('days'))
ts6[:] = stor_trans_BR_JS


