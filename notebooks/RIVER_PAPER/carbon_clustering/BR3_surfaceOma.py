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
def find_depth_deepalg(dp,prof,ts_y,ts_x, verbose = True):
    
    horizon = -9999
    water_depth = bathy['Bathymetry'][ts_y,ts_x]
    
    #case where no undersaturated water
    where_undersat = np.where(prof < 1); where_undersat = np.array(where_undersat)
    where_ssat = np.where((prof >= 1) & (prof < 1e19)); where_ssat = np.array(where_ssat)
    
    if(where_undersat.size == 0):
        if verbose:
            print('supersaturated to bottom')
        horizon = water_depth
        return horizon
    
    #case where only undersaturated water
    if(where_ssat.size == 0):
        if verbose:
            print('max aragonite found')
            print(np.max(prof[prof<1e19]))
        horizon = 0
        return horizon
    
    #case where last water cell is >1
    prof_water = prof[prof<1e20]
    last_water = prof_water[-1]
    if last_water >=1 :
        horizon = water_depth
        return horizon
    
    #case where there is some undersat water below all supersat water
    deepest_ssat = np.max(where_ssat)
    first_proper_undersat = (np.min(where_undersat[where_undersat>deepest_ssat]))
    horizon = (dp[first_proper_undersat]+dp[first_proper_undersat-1])/2
    
    return(horizon)

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
        toma = (tdat['model_output']['OmA'][0,ts_y,ts_x])
        #print(toma)
        daily_omahor[day] = toma #find_depth_deepalg(dp,toma,ts_y,ts_x, verbose = False)

    omah = xr.Dataset({'arag_0m':(['t'], daily_omahor)})
    stn_name = './ncs/' +  run_name + '_' + str(stn)  + 'asat_surf_sp' + str(spacing)+ '.nc'
    omah.to_netcdf(stn_name)
    

