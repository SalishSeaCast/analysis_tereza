yr = 2014
yrday = 365

import arrow
import pickle
import netCDF4 as nc
import numpy as np

tnam = f"TEMP_{yr}_JDF.pkl"
tnam2 = f"SAL_{yr}_JDF.pkl"

TEMP_ar = np.zeros([yrday,40,75])
SAL_ar = np.zeros([yrday,40,75])

start =f'{yr}-01-01'
end =f'{yr}-12-31'

start_run = arrow.get(start)
end_run = arrow.get(end)
arrow_array2 = []

for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array2.append(r)

for i in range(0,365):

    tdate = arrow_array2[i][0]
    yyyy = tdate.format('YYYY')
    yy = tdate.format('YY')
    mm = tdate.format('MM')
    dd = tdate.format('DD')
    mon = tdate.format('MMM')
    mon = mon.lower()
    verbstr = f'{dd}{mon}{yy}'
    print(verbstr)
    
    grid_nc = nc.Dataset(f'/results2/SalishSea/hindcast.201905/{dd}{mon}{yy}/SalishSea_1h_{yyyy}{mm}{dd}_{yyyy}{mm}{dd}_grid_T.nc')

    T = grid_nc['votemper'][:,:,361:(361+75),23]
    T_day = np.nanmean(T, axis = 0)
    TEMP_ar[i,:,:] = T_day
    
    S = grid_nc['vosaline'][:,:,361:(361+75),23]
    S_day = np.nanmean(S, axis = 0)
    SAL_ar[i,:,:] = S_day

pickle.dump(TEMP_ar, open(tnam, 'wb'))
pickle.dump(SAL_ar, open(tnam2, 'wb'))
