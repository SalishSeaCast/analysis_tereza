from numpy import *
from scipy import *

import netCDF4 as nc
import numpy as np
import scipy as sp
import cmocean
import glob
import seawater
import picklr2


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



NO3_mod = np.zeros_like(dic)
stn_check = np.zeros_like(dic)
js_f = np.zeros_like(dic)
is_f = np.zeros_like(dic)
ds_f = np.zeros_like(dic)
ws_f = np.zeros_like(dic)
js_ref = np.zeros_like(dic)
is_ref = np.zeros_like(dic)

datname = 'grid'
fname = 'vosaline'
firstmo = 1
yr = '18'


for i in range(0,len(dic)):
    if (i%10==0):
        print(i)
    t_NO3_mod, t_stn_with_nans, j_f, i_f, d_f, w_f, j_ref, i_ref = \
    picklr2.find_model_dat(i, fname, datname, firstmo, yr )
    
    NO3_mod[i] = t_NO3_mod
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

ncname = 'sal_RFACT_fy_2018.nc'
f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
g.createDimension('ref', len(NO3_mod))
ts = g.createVariable('mod_point','f4',('ref'))
ts[:] = NO3_mod[:]
ts2 = g.createVariable('mod_y','f4',('ref'))
ts2[:] = js_ref[:]
ts3 = g.createVariable('mod_x','f4',('ref'))
ts3[:] = is_ref[:]
f.close()