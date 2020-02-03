start = '2013-01-01'
end = '2013-12-31'
noday = 365

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import scipy as sp
import warnings
warnings.filterwarnings('ignore')
import datetime as dt
import glob
import arrow
import time
from datetime import datetime

def make_nclen(start,end):
    fn_ar_carp = []
    fn_ar_grid = []
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
        nc_sens_carp = '/results/SalishSea/hindcast.201905/' + ddmmmyy + '/*carp*.nc'
        tnc_sens_carp = glob.glob(nc_sens_carp)
        fn_ar_carp.append(tnc_sens_carp[0])
        nc_sens_grid = '/results/SalishSea/hindcast.201905/' + ddmmmyy + '/*grid_T*.nc'
        tnc_sens_grid = glob.glob(nc_sens_grid)
        fn_ar_grid.append(tnc_sens_grid[0])
        
    return fn_ar_carp, fn_ar_grid

def make_fname_ar(start,end,tstr):
    fn_ar = []
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    dayslen = len(arrow_array)
    for i in range(0,dayslen):
        tdate = arrow_array[i][0]
        ymd = tdate.format('YYYYMMDD')
        fname = tstr + ymd +'.nc'
        fn_ar.append(fname)
    
    return fn_ar

fn_ar_carp, fn_ar_grid = make_nclen(start,end)
tstr = 'grid_1d_'
avgdnc_ar_carp = make_fname_ar(start,end,tstr)

for i in range(0,noday):
    
    #should really be tgrid 
    tcarp = fn_ar_grid[i]
    fn = avgdnc_ar_carp[i]
    print(fn)
    
    t = time.time()
    carp = nc.Dataset(tcarp)
    SAL = carp['vosaline'][:]
    SAL_d = np.nanmean(SAL, axis = 0)
    TEMP = carp['votemper'][:]
    TEMP_d = np.nanmean(TEMP, axis = 0)
    t2 = time.time()
    print(t2-t)

    tdir = '/data/tjarniko/results/hindcast.201905_dayavg/'
    ncname = tdir + fn

    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
    g.createDimension('depths', 40)
    g.createDimension('ydir',898)
    g.createDimension('xdir',398)

    ts = g.createVariable('SAL','f4',('depths','ydir','xdir'))
    ts[:] = SAL_d
    ts = g.createVariable('TEMP','f4',('depths','ydir','xdir'))
    ts[:] = TEMP_d

    f.close()
