from __future__ import print_function
from numpy import *
from scipy import *
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

#ncname_BR = 'MASSBAL_BR_2015_2_fullyear.nc'
#sdir = 'BR_2nd_2015'
#start = '2015-01-01'
#end = '2015-12-31'
#st = dt.datetime(2015,1,1)
#en = dt.datetime(2015,12,31)
#create_massbal_nc(ncname_BR, sdir, start, end, st, en)

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
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        
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
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*' +'dian_U_' + ymd + '-' + ymd + '.nc'
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
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*' +'dian_V_' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        sens_ar.append(tnc_sens[0])

    return sens_ar

def calculate_total_C(files, size_domain):
    stor_mol = np.zeros(len(files))
    print('calculating total carbon')
    i = 0
    for f in files:
        if i%50 == 0:
            print(i)
        G = nc.Dataset(f)
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
    print('calculating surface carbon')
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

def calculate_surface20_C(files, size_domain_20):
    stor_mol = np.zeros(len(files))
    print('calculating surface 20 m carbon')
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

def calculate_20_100_C(files, size_domain_20_100):
    stor_mol = np.zeros(len(files))
    print('calculating carbon 20-100 m')
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

def calculate_100_deep_C(files, size_domain_deep):
    stor_mol = np.zeros(len(files))
    print('calculating carbon from 100m-deep')
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
    print('calculating air-sea flux')
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
    print('calculating transports JDF')
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
    print('calculating transports JS')
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


def create_massbal_nc(ncname, sdir, start, end, st, en, calc_JS_trans):
    
    if calc_JS_trans == 1:
        print('calculating JS transports')
    if calc_JS_trans != 1:
        print('not calculating JS transports')

    y_st = st.timetuple().tm_yday
    print(y_st)
    y_en = en.timetuple().tm_yday
    print(y_en)
    ts_BR = np.arange(y_st,y_en+1,1)
    days_in = np.size(ts_BR)

    BR_ar = make_nclen(start,end,'carp', sdir)
    BR_ar_tp = make_nclen_transport(start,end, sdir)
    if calc_JS_trans == 1:
        BR_ar_tpV = make_nclen_transport_V(start,end, sdir)

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
    
    stor_mol = calculate_total_C(BR_ar, size_domain)
    stor_mol_surf = calculate_surface_C(BR_ar, surfa)
    stor_mol_20 = calculate_surface20_C(BR_ar, size_domain_20)
    stor_mol_20_100 = calculate_20_100_C(BR_ar, size_domain_20_100)
    stor_mol_deep = calculate_100_deep_C(BR_ar, size_domain_deep)
    stor_flx = calculate_flux(BR_ar, surfa)
    stor_trans_JDF = calculate_transports(BR_ar_tp)
    if calc_JS_trans == 1:
        stor_trans_JS = calculate_transports_JS(BR_ar_tpV)


    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
    g.createDimension('days', days_in)
    ts = g.createVariable('stor_mol','f4',('days'))
    ts[:] = stor_mol
    ts2 = g.createVariable('stor_mol_surf','f4',('days'))
    ts2[:] = stor_mol_surf
    ts3 = g.createVariable('stor_mol_20','f4',('days'))
    ts3[:] = stor_mol_20
    ts3b = g.createVariable('stor_mol_20_100','f4',('days'))
    ts3b[:] = stor_mol_20_100
    ts3c = g.createVariable('stor_mol_deep','f4',('days'))
    ts3c[:] = stor_mol_deep
    ts4 = g.createVariable('stor_flx','f4',('days'))
    ts4[:] = stor_flx
    ts5 = g.createVariable('stor_trans_JDF','f4',('days'))
    ts5[:] = stor_trans_JDF
    if calc_JS_trans == 1:
        ts6 = g.createVariable('stor_trans_JS','f4',('days'))
        ts6[:] = stor_trans_JS
        
    f.close()



