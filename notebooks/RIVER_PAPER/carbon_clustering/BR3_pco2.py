import sys
sys.path.append('./extraction_scripts/')
import arrow
import xarray as xr
import numpy as np
import glob
import netCDF4 as nc
import map_fxn as mf

bathy = nc.Dataset('/data/tjarniko/MEOPAR/grid/bathymetry_201702.nc')
############oma eq###############


str_to_oma = '/data/tjarniko/results/BASERUN_EXP/Oma_calc/BR3_OmA_plus_'
run_name = 'BR3'
stn_b = 0
stn_e = 580

start ='2015-01-01'
end ='2015-12-31'

start_run = arrow.get(start)
end_run = arrow.get(end)

arrow_array1 = []

for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array1.append(r)
dayslen = len(arrow_array1)

dayar = []
for i in range(0,dayslen):

    tdate = arrow_array1[i][0]
    ymd = tdate.format('YYYYMMDD')
    tstr = glob.glob(str_to_oma +ymd+'*.nc')
    dayar.append(tstr[0])    
    
bath = '/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc'
grid = mf.import_bathy(bath)
fmask = (grid.fmask[0,0,:,:])
spacing = 10
stn_x, stn_y = mf.make_stns(spacing)
d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)    


tdp = grid['gdept_0'][:]
dp = np.squeeze(tdp)
## start extraction per station
for stn in range(stn_b,stn_e):

    print('station is: ' ,str(stn))
    print('x is :', d_stn_x[stn])
    print('y is :', d_stn_y[stn])

    ts_x = d_stn_x[stn]
    ts_y = d_stn_y[stn]
    
    daily_omahor = np.zeros(len(dayar))

    for day in range(0,len(dayar)):
        if day%20 == 0:
            print(day)
        tdat = nc.Dataset(dayar[day])
        daily_omahor[day] = (tdat['model_output']['pCO2'][ts_y,ts_x])
        #print(toma)

        omah = xr.Dataset({'pCO2':(['t'], daily_omahor)})
        stn_name = './ncs/' +  run_name + '_' + str(stn)  + 'surfpco2_sp' + str(spacing)+ '.nc'
        omah.to_netcdf(stn_name)
    

