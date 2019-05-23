from __future__ import print_function
from numpy import *
from scipy import *
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import scipy as sp
import cmocean
import glob
import seawater
""
from salishsea_tools import (
    nc_tools,
    viz_tools,
    geo_tools,
    tidetools
)

infil = loadtxt('../DATASETS/grl2016_edit2.txt')
infil_cor = loadtxt('../DATASETS/grl2016_nu.txt')

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

def find_model_dat(w, fname, datname ):
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
    q2 = nc.Dataset('/results2/SalishSea/nowcast-green.201806/01jan18/SalishSea_1h_20180101_20180101_grid_T.nc')
    dep=q2.variables['deptht']
    d_mod = dep[:]  
    
    #months, to make a path to model data
    mons = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    month = mon[w]

    #TJ check for month in may or later, findable station, and good quality flag, and dic>0 (should be redundant)
    if ((month>=5) & (~np.isnan(j)) & (~np.isnan(i)) & ((dic_qf[w]==2)|(dic_qf[w]==6)) &(dic[w]>0)):
        
        #from month and day, make string to 
        mon_name = str(mons[int(month)-1])
        day_name = str(int(day[w]))
        if len(day_name) ==1:
            day_name = '0'+day_name
        modpath = day_name+mon_name+'17/'

        bigpath = '/results2/SalishSea/hindcast.201812_annex/'
        if datname == 'carp':
            ncpath = 'SalishSea_1d_*_carp_T.nc'
        if datname == 'ptrc':
            ncpath = 'SalishSea_1d_*_ptrc_T.nc'
        else:
            ncpath = 'SalishSea_1d_*_grid_T.nc'
        q = glob.glob(bigpath+modpath+ncpath)
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
        

sal_mod = np.zeros_like(dic)
stn_check = np.zeros_like(dic)
js_f = np.zeros_like(dic)
is_f = np.zeros_like(dic)
ds_f = np.zeros_like(dic)
ws_f = np.zeros_like(dic)
js_ref = np.zeros_like(dic)
is_ref = np.zeros_like(dic)

datname = 'grid'
fname = 'votemper'
for i in range(0,len(dic)):
    if (i%10==0):
        print(i)
    t_sal_mod, t_stn_with_nans, j_f, i_f, d_f, w_f, j_ref, i_ref = find_model_dat(i, fname, datname )
    sal_mod[i] = t_sal_mod
    stn_check[i] = t_stn_with_nans
    js_f[i]=j_f
    is_f[i]=i_f
    ds_f[i]=d_f
    ws_f[i]=w_f
    js_ref[i]=j_ref
    is_ref[i]=i_ref


itt = is_f[is_f<9999]
jtt = js_f[js_f<9999]
dtt = ds_f[ds_f<9999]
wtt = ws_f[ws_f<9999]
print(wtt)
print('**')
itt2 = is_f[is_f<9999]
jtt2 = js_f[is_f<9999]
dtt2 = ds_f[is_f<9999]
wtt2 = ws_f[is_f<9999]
print(wtt2)


#t_dic_mod, t_stn_with_nans, j_f, i_f, d_f, w_f, j_ref, i_ref

ncname = 'temp_may012017_end2017.nc'
f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
g.createDimension('ref', len(sal_mod))
ts = g.createVariable('mod_point','f4',('ref'))
ts[:] = sal_mod[:]
ts2 = g.createVariable('mod_y','f4',('ref'))
ts2[:] = js_ref[:]
ts3 = g.createVariable('mod_x','f4',('ref'))
ts3[:] = is_ref[:]
f.close()