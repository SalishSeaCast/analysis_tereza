import arrow
import pickle
import netCDF4 as nc
import numpy as np

yr = 2016
yrday = 366
trans_sal_2015_mol_per_gridcell = np.zeros([yrday,40,75])
tnam = f"trans_sal_{yr}_mol_per_gridcell.pkl"

start =f'{yr}-01-01'
end =f'{yr}-12-31'

start_run = arrow.get(start)
end_run = arrow.get(end)
arrow_array2 = []

for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array2.append(r)

for i in range(0,yrday):

    tdate = arrow_array2[i][0]
    yyyy = tdate.format('YYYY')
    yy = tdate.format('YY')
    mm = tdate.format('MM')
    dd = tdate.format('DD')
    mon = tdate.format('MMM')
    mon = mon.lower()
    #print(mon)
    verbstr = f'{dd}{mon}{yy}'
    print(verbstr)
    
    uflx_nc = nc.Dataset(f'/results2/SalishSea/hindcast.201905/{dd}{mon}{yy}/SalishSea_1h_{yyyy}{mm}{dd}_{yyyy}{mm}{dd}_flux_U.nc')
    grid_nc = nc.Dataset(f'/results2/SalishSea/hindcast.201905/{dd}{mon}{yy}/SalishSea_1h_{yyyy}{mm}{dd}_{yyyy}{mm}{dd}_grid_T.nc')
   

    # transport in kg/s at the i = 23, j = 361:361+75 boundary
    u_transport = uflx_nc['u_masstr'][:,:,:,0]
    # sal concentration there
    sal = grid_nc['vosaline'][:,:,361:(361+75),23]

    #calculate sal in kg/hour
    #kg/g* kg/s * s/ hr to get flux per hour
    C_hr =  1e-3 * sal * u_transport * (60*60)
    #sum hourly flux to get umol C/day
    C_day = np.nansum(C_hr, axis = 0)
    trans_sal_2015_mol_per_gridcell[i,:,:] = C_day
    

pickle.dump(trans_sal_2015_mol_per_gridcell, open(tnam, 'wb'))