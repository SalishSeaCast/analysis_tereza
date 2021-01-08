start = '2016-03-01'
end = '2016-08-31'
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
    fn_ar_ptrc = []
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
        nc_sens_ptrc = '/results2/SalishSea/hindcast.201905/' + ddmmmyy + '/*ptrc*.nc'
        tnc_sens_ptrc = glob.glob(nc_sens_ptrc)
        fn_ar_ptrc.append(tnc_sens_ptrc[0])
        nc_sens_grid = '/results2/SalishSea/hindcast.201905/' + ddmmmyy + '/*grid*.nc'
        tnc_sens_grid = glob.glob(nc_sens_grid)
        fn_ar_grid.append(tnc_sens_grid[0])

    return fn_ar_ptrc, fn_ar_grid

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

fn_ar_ptrc, fn_ar_grid = make_nclen(start,end)
tstr = 'ptrc_phyto_1d_'
avgdnc_ar_ptrc = make_fname_ar(start,end,tstr)

for i in range(0,noday):

    tptrc = fn_ar_ptrc[i]
    fn = avgdnc_ar_ptrc[i]
    print(fn)

    t = time.time()
    ptrc = nc.Dataset(tptrc)
#     print(ptrc)
    diat = ptrc['diatoms'][:]
    diat_d = np.nanmean(diat, axis = 0)
    flag = ptrc['flagellates'][:]
    flag_d = np.nanmean(flag, axis = 0)
    cili = ptrc['ciliates'][:]
    cili_d = np.nanmean(cili, axis = 0)
    ptrc.close()
    t2 = time.time()
    print(t2-t)

    tdir = '/data/tjarniko/results/hindcast.201905_dayavg_phyto/'
    ncname = tdir + fn

    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
    g.createDimension('depths', 40)
    g.createDimension('ydir',898)
    g.createDimension('xdir',398)

    ts = g.createVariable('diatoms','f4',('depths','ydir','xdir'))
    ts[:] = diat_d
    ts = g.createVariable('flagellates','f4',('depths','ydir','xdir'))
    ts[:] = flag_d
    ts = g.createVariable('ciliates','f4',('depths','ydir','xdir'))
    ts[:] = cili_d

    f.close()
