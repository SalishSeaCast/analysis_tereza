import sys
sys.path.append('./KEY_PAPERFIGURES/extraction_scripts/')
import arrow
import xarray as xr
import numpy as np
import glob
import netCDF4 as nc
import map_fxn as mf

bathy = nc.Dataset('/data/tjarniko/MEOPAR/grid/bathymetry_201702.nc')
############oma eq###############


str_to_W = '/data/tjarniko/avg/*'
run_name = 'BR3'
stn_b = 300
stn_e = 400

start ='2013-01-01'
end ='2013-12-31'

start_run = arrow.get(start)
end_run = arrow.get(end)

arrow_array1 = []

for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array1.append(r)
dayslen = len(arrow_array1)

gridW = []
for i in range(0,dayslen):

    tdate = arrow_array1[i][0]
    ymd = tdate.format('YYYYMMDD')
    tstr = glob.glob(str_to_W +ymd+'*grid_W.nc')
    #print(tstr)
    gridW.append(tstr[0])    

carpT = []
str_to_C = '/results2/SalishSea/hindcast.201905/*/*'
for i in range(0,dayslen):

    tdate = arrow_array1[i][0]
    ymd = tdate.format('YYYYMMDD')
    tstr = (str_to_C+ymd+'*carp_T.nc')#'*carp_T.nc'
    #print(tstr)
    tstr = glob.glob(str_to_C +ymd+'*carp_T.nc')
    #print(tstr)
    #print(tstr)
    carpT.append(tstr[0])       
    
bath = '/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc'
grid = mf.import_bathy(bath)
fmask = (grid.fmask[0,0,:,:])
spacing = 10
stn_x, stn_y = mf.make_stns(spacing)
d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)  
print('complete')

year = '2013'
LL=0.1 

def getVEDEuph(ii,jj,LL,fileK,filePAR):
    #print(filePAR)
#     print(ii)
#     print(jj)
    dailyPAR=np.mean(filePAR.variables['PAR'][:,:,jj,ii],0)
    kk=np.sum(dailyPAR>LL*dailyPAR[0])
    #print(fileK)
    kprof=(fileK.variables['vert_eddy_diff'][:,:,jj,ii])
    #print(kprof)
    new_kprof = kprof[0,:]
    return new_kprof[kk]


for stn in range(stn_b,stn_e):

    print('station is: ' ,str(stn))
    print('x is :', d_stn_x[stn])
    print('y is :', d_stn_y[stn])

    ts_x = d_stn_x[stn]
    ts_y = d_stn_y[stn]

    daily_omahor = np.zeros(dayslen)

    for day in range(0,dayslen):

        carp = nc.Dataset(carpT[day])
        W = nc.Dataset(gridW[day])
        #print(gridW[day])
        vedeuph = getVEDEuph(ts_x,ts_y,LL,W,carp)
        if day%5 == 0:
            print(day)
            print(vedeuph)
        daily_omahor[day] = vedeuph
        #print(toma)

        ved = xr.Dataset({'daily_ved':(['t'], daily_omahor)})
        stn_name = '/data/tjarniko/MEOPAR/analysis_tereza/notebooks/CLUSTER_PAPER/CLEAN/NC_HINDCAST/'\
        + str(year) + '/VED_TS/stn_' + str(stn)  + 'FOTOVED_sp' + str(spacing)+ '.nc'
        ved.to_netcdf(stn_name)
        