import arrow
import pickle
import netCDF4 as nc
import numpy as np

yr = 2019
yrday = 365
trans_DIC_2019_mol_per_gridcell = np.zeros([365,40,75])
anom_DIC = 2050
tnam = f"trans_anomDIC_{anom_DIC}_2019_mol_per_gridcell.pkl"

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
    #print(mon)
    verbstr = f'{dd}{mon}{yy}'
    print(verbstr)
    
    uflx_nc = nc.Dataset(f'/results2/SalishSea/hindcast.201905/{dd}{mon}{yy}/SalishSea_1h_{yyyy}{mm}{dd}_{yyyy}{mm}{dd}_flux_U.nc')
    carp_nc = nc.Dataset(f'/results2/SalishSea/hindcast.201905/{dd}{mon}{yy}/SalishSea_1h_{yyyy}{mm}{dd}_{yyyy}{mm}{dd}_carp_T.nc')
   

    # transport in kg/s at the i = 23, j = 361:361+75 boundary
    u_transport = uflx_nc['u_masstr'][:,:,:,0]
    # DIC concentration there
    DIC = carp_nc['dissolved_inorganic_carbon'][:,:,361:(361+75),23]- anom_DIC
    #potential density sigma theta there in kg/L
    sigma0 = (carp_nc['sigma_theta'][:,:,361:(361+75),23]+1000)/1000
    #inverse of pot. dense in L/kg
    inv_sigma0 = 1/sigma0

    #calculate C mol per hour
    #mol/umol * umol/L * L/kg * kg/s * s/ hr to get flux per hour
    C_hr =  1e-6 * DIC * inv_sigma0 * u_transport * (60*60)
    #sum hourly flux to get umol C/day
    C_day = np.nansum(C_hr, axis = 0)
    trans_DIC_2019_mol_per_gridcell[i,:,:] = C_day
    

pickle.dump(trans_DIC_2019_mol_per_gridcell, open(tnam, 'wb'))