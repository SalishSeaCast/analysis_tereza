from numpy import *
from scipy import *

import netCDF4 as nc
import numpy as np
import scipy as sp
import cmocean
import glob
import seawater
import picklr_forhindcast201905
import time

from salishsea_tools import (
    nc_tools,
    viz_tools,
    geo_tools,
    tidetools)

start = time.time()
print(start)

#constants (path and yrs)
bigpath = '/data/tjarniko/results/hindcast.201905_dayavg_DIC-TA-T-S/'
firstmos = [1,1,1,1,1]
yrs = ['13','14','15','16','17']

#uncomment whichever variable you are creating a dataset for
#TA
# ncnames = ['TA_hindcast201905_GRL_comp_2013.nc','TA_hindcast201905_GRL_comp_2014.nc',
#           'TA_hindcast201905_GRL_comp_2015.nc','TA_hindcast201905_GRL_comp_2016.nc',
#           'TA_hindcast201905_GRL_comp_2017.nc']        
# datnames = ['carp','carp','carp','carp','carp']
# fnames = ['TA','TA','TA','TA','TA']

#DIC
# ncnames = ['DIC_hindcast201905_GRL_comp_2013.nc','DIC_hindcast201905_GRL_comp_2014.nc',
#           'DIC_hindcast201905_GRL_comp_2015.nc','DIC_hindcast201905_GRL_comp_2016.nc',
#           'DIC_hindcast201905_GRL_comp_2017.nc']        
# datnames = ['carp','carp','carp','carp','carp']
# fnames = ['DIC','DIC','DIC','DIC','DIC']

#sal
# ncnames = ['sal_hindcast201905_GRL_comp_2013.nc','sal_hindcast201905_GRL_comp_2014.nc',
#           'sal_hindcast201905_GRL_comp_2015.nc','sal_hindcast201905_GRL_comp_2016.nc',
#           'sal_hindcast201905_GRL_comp_2017.nc']        
# datnames = ['grid','grid','grid','grid','grid']
# fnames = ['SAL','SAL','SAL','SAL','SAL']

# #temp
ncnames = ['temp_hindcast201905_GRL_comp_2013.nc','temp_hindcast201905_GRL_comp_2014.nc',
          'temp_hindcast201905_GRL_comp_2015.nc','temp_hindcast201905_GRL_comp_2016.nc',
          'temp_hindcast201905_GRL_comp_2017.nc']        
datnames = ['grid','grid','grid','grid','grid']
fnames = ['TEMP','TEMP','TEMP','TEMP','TEMP']


#load grl data for comparison
infil = loadtxt('../Datasets/grl2016_nu.txt')

crid= infil[:,0]
ln = infil[:,2]
stn = infil[:,3]
mon = infil[:,4]
day = infil[:,5]
lat_or = infil[:,6]
lon_or = infil[:,7]
lat = infil[:,6]
lon = infil[:,7]
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



for t_var in range(0,len(datnames)):
    
    ncname = ncnames[t_var]
    datname = datnames[t_var] 
    fname = fnames[t_var]
    firstmo = firstmos[t_var]
    yr = yrs[t_var]
    
    print(ncname)
    
    stn_check = np.zeros_like(dic)
    js_f = np.zeros_like(dic)
    is_f = np.zeros_like(dic)
    ds_f = np.zeros_like(dic)
    ws_f = np.zeros_like(dic)
    js_ref = np.zeros_like(dic)
    is_ref = np.zeros_like(dic)
    var_mod = np.zeros_like(dic)

    for i in range(0,len(dic)):
        if (i%10==0):
            print(i)
        t_var_mod, t_stn_with_nans, j_f, i_f, d_f, w_f, j_ref, i_ref = \
        picklr_forhindcast201905.find_model_dat(i, fname, datname, firstmo, yr, bigpath )

        var_mod[i] = t_var_mod
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



    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
    g.createDimension('ref', len(var_mod))
    ts = g.createVariable('mod_point','f4',('ref'))
    ts[:] = var_mod[:]
    ts2 = g.createVariable('mod_y','f4',('ref'))
    ts2[:] = js_ref[:]
    ts3 = g.createVariable('mod_x','f4',('ref'))
    ts3[:] = is_ref[:]
    f.close()

end = time.time()
print('time of code run')
print(end-start)