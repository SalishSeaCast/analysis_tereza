from numpy import *
from scipy import *

import netCDF4 as nc
import numpy as np
import scipy as sp
import cmocean
import glob
import seawater
import arrow



from salishsea_tools import (
    nc_tools,
    viz_tools,
    geo_tools,
    tidetools
)

infil = loadtxt('/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/DATASETS/grl2016_edit2.txt')
infil_cor = loadtxt('/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/DATASETS/grl2016_nu.txt')

crid= infil[:,0]
ln = infil[:,2]
stn = infil[:,3]
mon = infil[:,4]
day = infil[:,5]
lat_or = infil[:,6]
lon_or = infil[:,7]
lat = infil_cor[:,6]
lon = infil_cor[:,7]
P = infil[:,8]
T = infil[:,9]
S = infil[:,10]
ox = infil[:,11]
ox_qf = infil[:,12]
dic = infil[:,13]
alk = infil[:,15]
dic_qf = infil[:,14]
alk_qf = infil[:,16]
no3 = infil[:,17]
no3_qf = infil[:,18]
si = infil[:,19]
si_qf = infil[:, 20]
po4 = infil[:,21]
po4_qf = infil[:, 22]

mesh = nc.Dataset('/data/tjarniko/MEOPAR/grid/bathymetry_201702.nc')
nav_lon = mesh.variables['nav_lon'][:]
nav_lat = mesh.variables['nav_lat'][:]
bathy = mesh.variables['Bathymetry'][:]


def find_model_dat(w, fname, datname, firstmo, yr, path_to_ncs):
    "Usage: "
    model_dat = np.nan
    
    #if the stn has lat, lons not findable in model (ie, out of model domain) this will change to the # of the station
    stn_with_nans = 9999
    
    #if the observations have data but the model has a 0, these (_f = faulty) record model lat,lon,depth
    w_f = 9999
    d_f = 9999
    i_f = 9999
    j_f = 9999
    
    #find closest model point, check that it's not nan, if it is, record it
    j, i = geo_tools.find_closest_model_point(lon[w],lat[w],nav_lon,nav_lat)
    j_ref = j
    i_ref = i
    if (np.isnan(i) | np.isnan(j)):
        stn_with_nans = stn[w]
    
    #extract model depths from model, for comparison with obs
    q2 = nc.Dataset('/results2/SalishSea/hindcast.201905/07apr10/SalishSea_1h_20100407_20100407_grid_T.nc')
    dep=q2.variables['deptht']
    d_mod = dep[:]  

    month = mon[w]

    #TJ check for month in first month or later, findable station, and good quality flag, and dic>0 (should be redundant)
    if ((month>=firstmo) & (~np.isnan(j)) & (~np.isnan(i)) & ((dic_qf[w]==2)|(dic_qf[w]==6)) &(dic[w]>0)):
        
        #from month and day, make string to 
        #mon_name = str(mons[int(month)-1])
        mon_name = str(int(mon[w]))
        if len(mon_name) ==1:
            mon_name = '0'+mon_name
        day_name = str(int(day[w]))
        if len(day_name) ==1:
            day_name = '0'+day_name
        #modpath = day_name+mon_name+ yr+'/'
        start = '20' + yr + '-' +mon_name+ '-' +day_name 
        tdate =arrow.get(start)
        ymd = tdate.format('YYYYMMDD')
        print(ymd)

        bigpath = path_to_ncs
        #print(datname)
        ncpath = '/SKOG_1d_*'+ datname +'_T_' + ymd + '-' + ymd + '.nc'
        q = glob.glob(bigpath+ncpath)
        
        #depths of observations
        dtt = P[w]
        
        #find closest depth in model - print if necesary
        close = np.abs(d_mod-dtt)
        t_depth = np.argmin(close)
        depth_mod = d_mod[t_depth]
#         print('DEPTH IND: '+str(t_depth))
#         print('DEPTH : '+str(depth_mod))
        
        #open dataset, record relevant model DIC
        qnc = nc.Dataset(q[0])
        #print(qnc)
        #print(fname)
        model_dat = qnc.variables[fname][0,t_depth,j,i]
        # flag spots where model gives 0s. 
        if (model_dat == 0):
        # we have     
#             print('*******')
#             print('index of grl data: '+str(w))
#             print('depth in grl data: '+str(P[w]))
#             print('j index: '+str(j))
#             print('i index: '+str(i))
#             print('DEPTH IND: '+str(t_depth))
#             print('DEPTH : '+str(depth_mod))
            w_f = w
            j_f = j
            i_f = i
            d_f = t_depth
            
#         print('Data dic: '+str( dic[w]))
#         print('Data dic qf: '+str( dic_qf[w]))
#         print('Model dic '+str(model_dat))
    
    return model_dat, stn_with_nans, j_f, i_f, d_f, w_f, j_ref, i_ref