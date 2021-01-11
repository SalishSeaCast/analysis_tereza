import arrow
import netCDF4 as nc
import glob
import numpy as np
import time
import xarray as xr
import matplotlib.pyplot as plt
from salishsea_tools import visualisations as vis
import cmocean as cm

spr_st = 59; spr_e = 151; sum_st = 151; sum_e = 243;
#'2015-03-01' '2015-06-01'
#'2015-06-01' '2015-08-31

#/results2/SalishSea/hindcast.201905/'
yr = '2013'
start =f'{yr}-06-01'
end =f'{yr}-08-31'

start_run = arrow.get(start)
end_run = arrow.get(end)

arrow_array = []
for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array.append(r)
dayslen = len(arrow_array)
print(dayslen)
    
diat = np.zeros([dayslen,40,898,398])
flag = np.zeros([dayslen,40,898,398])
cili = np.zeros([dayslen,40,898,398])

for i in range(0,dayslen):
    

    tdate = arrow_array[i][0]
    
    mm = tdate.format('MM')
    dd = tdate.format('DD')
    yy = tdate.format('YYYY')
    d_str = f'{yy}{mm}{dd}'
    ptrc_str = f'/data/tjarniko/results/hindcast.201905_dayavg_phyto/ptrc_phyto_1d_{d_str}.nc'
    t_dat = glob.glob(ptrc_str)
    print(t_dat[0])
    t0 = time.time()
    w = nc.Dataset(t_dat[0])
    diat[i,:,:,:] = w['model_output']['diatoms'][:]
    flag[i,:,:,:] = w['model_output']['flagellates'][:]
    cili[i,:,:,:] = w['model_output']['ciliates'][:]
    w.close()
    t1 = time.time()
    total = t1-t0
    #print(f'seconds taken: {str(total)}')
    
cili_s = np.nanmean(cili, axis = 0)
flag_s = np.nanmean(flag, axis = 0)
diat_s = np.nanmean(diat, axis = 0)

ncname = './pkls/2013_SUM_phyto.nc'


f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('model_output')
g.createDimension('depths', 40)
g.createDimension('ydir',898)
g.createDimension('xdir',398)

ts = g.createVariable('diatoms','f4',('depths','ydir','xdir'))
ts[:] = diat_s
ts = g.createVariable('flagellates','f4',('depths','ydir','xdir'))
ts[:] = flag_s
ts = g.createVariable('ciliates','f4',('depths','ydir','xdir'))
ts[:] = cili_s

f.close()
