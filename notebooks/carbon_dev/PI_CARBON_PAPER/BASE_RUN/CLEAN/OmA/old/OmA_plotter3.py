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


start1 = '2015-01-01'
end1 = '2015-06-30'
sdir_preind = '/data/tjarniko/results/BASERUN_EXP/PI_3rd_2015/ncs/'
sdir_br = '/data/tjarniko/results/BASERUN_EXP/BR_2nd_2015/ncs/'
figstr = 'Jan1_Jun30_OmA.png'
ncname = 'Jan1_Jun30_OmA.nc'


def make_nclen(start,end,ftype, sdir):
    date_ar = []
    sens_ar = []
    doy_ar = []
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
        dddd = tdate.format('DDDD')
        nc_sens = sdir + '/SKOG_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        sens_ar.append(tnc_sens[0])
        date_ar.append(ddmmmyy)
        doy_ar.append(dddd)
    return date_ar, sens_ar, doy_ar

def oned_moxy(tsal, ttemp, tdic, tta, pres_atm, depth_this):

    size_box = np.shape(tdic)
    size_0 = size_box[0]
    size_1= size_box[1]
    size_2 = size_box[2]


    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    #tdepth = np.zeros_like(tsra)
    tpressure[:] = pres_atm
    tdepth = np.ravel(depth_this)
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504
    ttera_is = gsw.t_from_CT(tsra,ttera,tzero)
    print('beginning mocsy')
    response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tdepth, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup
    print('finished mocsy')

    pHr = pH.reshape(size_0,size_1,size_2)
    OmAr = OmegaA.reshape(size_0,size_1,size_2)
    pco2r = pco2.reshape(size_0,size_1,size_2)
    
    return pHr, OmAr, pco2r

def make_avgthwg_plot_OmA(start,end,sdir_PI, sdir_BR, figstr):

    yr_s = int(start[0:4])
    mon_s = int(start[5:7])
    day_s = int(start[8:10])

    yr_e = int(end[0:4])
    mon_e = int(end[5:7])
    day_e = int(end[8:10])

    st = dt.datetime(yr_s,mon_s,day_s)
    en = dt.datetime(yr_e,mon_e,day_e)
    
    y_st = st.timetuple().tm_yday
    print(y_st)
    y_en = en.timetuple().tm_yday
    print(y_en)
    ts_BR = np.arange(y_st,y_en+1,1)
    days_in = len(ts_BR)

    dates_preind_carp, files_preind_carp, doy_preind = make_nclen(start,end,'carp', sdir_PI)
    dates_br_carp, files_br_carp, doy_br = make_nclen(start,end,'carp', sdir_BR)
    dates_preind_grid, files_preind_grid, doy_preind = make_nclen(start,end,'grid_T', sdir_PI)
    dates_br_grid, files_br_grid, doy_br = make_nclen(start,end,'grid_T', sdir_BR)
    
    dfile = nc.Dataset('/data/tjarniko/results/BASERUN_EXP/PI_3rd_2015/ncs/SKOG_1d_20150101_20150301_carp_T_20150101-20150101.nc')
    depths = dfile['deptht'][:]
    depth_broad = np.zeros([1,40,898,398])
    depth_broad2 = np.zeros([1,898,398])

    #expand_dims
    for i in range(0,40):
        depth_broad2[:] = depths[i]
        depth_broad[:,i,:,:] = depth_broad2

    mon3_dic_BR = np.zeros([days_in,40,898,398])
    mon3_dic_PI = np.zeros([days_in,40,898,398])
    mon3_ta_BR = np.zeros([days_in,40,898,398])
    mon3_ta_PI = np.zeros([days_in,40,898,398])
    mon3_temp_BR = np.zeros([days_in,40,898,398])
    mon3_temp_PI = np.zeros([days_in,40,898,398])
    mon3_sal_BR = np.zeros([days_in,40,898,398])
    mon3_sal_PI = np.zeros([days_in,40,898,398])
    mon3_OmA_BR = np.zeros([days_in,40,898,398])
    mon3_OmA_PI = np.zeros([days_in,40,898,398]) 
    
    for i in range (0,days_in):
        if i%5 ==0:
            print(i)
        test_br_carp = nc.Dataset(files_br_carp[i])
        test_pi_carp = nc.Dataset(files_preind_carp[i])
        test_br_grid = nc.Dataset(files_br_grid[i])
        test_pi_grid = nc.Dataset(files_preind_grid[i])
        t_dic_br = np.squeeze(test_br_carp['dissolved_inorganic_carbon'][:])
        t_dic_pi = np.squeeze(test_pi_carp['dissolved_inorganic_carbon'][:])
        t_ta_br = np.squeeze(test_br_carp['total_alkalinity'][:])
        t_ta_pi = np.squeeze(test_pi_carp['total_alkalinity'][:])
        t_sal_br = np.squeeze(test_br_grid['vosaline'][:])
        t_sal_pi = np.squeeze(test_pi_grid['vosaline'][:])        
        t_temp_br = np.squeeze(test_br_grid['votemper'][:])
        t_temp_pi = np.squeeze(test_pi_grid['votemper'][:])  
        
        mon3_dic_BR[i,:,:,:] = t_dic_br
        mon3_dic_PI[i,:,:,:] = t_dic_pi
        mon3_ta_BR[i,:,:,:] = t_ta_br
        mon3_ta_PI[i,:,:,:] = t_ta_pi
        mon3_sal_BR[i,:,:,:] = t_sal_br
        mon3_sal_PI[i,:,:,:] = t_sal_pi
        mon3_temp_BR[i,:,:,:] = t_temp_br
        mon3_temp_PI[i,:,:,:] = t_temp_pi

        pHr_pi, OmAr_pi, pco2r_pi = oned_moxy(t_sal_pi, t_temp_pi, t_dic_pi, t_ta_pi, 1, depth_broad)
        pHr_br, OmAr_br, pco2r_br = oned_moxy(t_sal_br, t_temp_br, t_dic_br, t_ta_br, 1, depth_broad)
        
        mon3_OmA_BR[i,:,:,:] = OmAr_br
        mon3_OmA_PI[i,:,:,:] = OmAr_pi
        
    mon3_OmA_BR_m = np.ma.masked_where(mon3_OmA_BR >= 1e10, mon3_OmA_BR)
    mon3_OmA_PI_m = np.ma.masked_where(mon3_OmA_PI >= 1e10, mon3_OmA_PI)
    OmAr_pi_av = np.mean(mon3_OmA_PI_m,axis = 0)
    OmAr_br_av = np.mean(mon3_OmA_BR_m,axis = 0)
    

    
    return mon3_OmA_BR_m, mon3_OmA_PI_m, days_in

mon3_OmA_BR_m, mon3_OmA_PI_m, days_in = make_avgthwg_plot_OmA(start1,end1,sdir_preind, sdir_br, figstr)

f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
#g.createDimension('days', len(NO3_mod))
g.createDimension('days', days_in)
g.createDimension('depths', 40)
g.createDimension('ys', 898)
g.createDimension('xs', 398)
ts = g.createVariable('OmAr_pi','f4',('days','depths','ys','xs'))
ts[:] = mon3_OmA_PI_m
ts2 = g.createVariable('OmAr_br','f4',('days','depths','ys','xs'))
ts2[:] = mon3_OmA_BR_m
f.close()

elapsed = timeit.default_timer() - start_time
print(elapsed)