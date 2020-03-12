start = '2015-01-01'
end = '2015-12-31'
## too much import things
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import scipy as sp
import warnings
warnings.filterwarnings('ignore')
import datetime as dt
import glob
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
import arrow
import gsw
from datetime import datetime

### 

###fxn definitions!

def make_nclen(start,end,ftype):
    fn_ar = []
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
        nc_sens = '/data/tjarniko/results/hindcast.201905_dayavg_DIC-TA-T-S/'+ ftype +'_1d_' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        fn_ar.append(tnc_sens[0])

    return fn_ar

def make_fname_ar(start,end):
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
        fname = 'OmA_plus_' + ymd +'.nc'
        fn_ar.append(fname)
    
    return fn_ar

def OmA_3D(grid,carp):
    tsal = grid['model_output']['SAL'][:,:,:]
    ttemp = grid['model_output']['TEMP'][:,:,:]
    tdic = carp['model_output']['DIC'][:,:,:]
    tta = carp['model_output']['TA'][:,:,:]
    
    test_LO = nc.Dataset('/results/forcing/LiveOcean/boundary_conditions/LiveOcean_v201905_y2018m01d01.nc')
    zlevels = (test_LO['deptht'][:])

    depths = np.zeros([40,898,398])

    for j in range(0,898):
        for i in range(0,398):
            depths[:,j,i] = zlevels
            
    tdepths = np.ravel(depths)
    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    tpressure[:] =1
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504
    ttera_is = gsw.t_from_CT(tsra,ttera,tzero)

    response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tdepths, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    pHr = pH.reshape(40,898,398)
    OmAr = OmegaA.reshape(40,898,398)
    OmCr = OmegaC.reshape(40,898,398)
    pco2r = pco2.reshape(40,898,398)
    
    return pHr, OmAr, OmCr, pco2r

############let's get extractin bro

carp_ar = make_nclen(start,end,'carp')
grid_ar = make_nclen(start,end,'grid')
fn_ar = make_fname_ar(start,end)

for i in range(0,365):
    
    tcarp = carp_ar[i]
    print(tcarp)
    carp = nc.Dataset(tcarp)
    tgrid = grid_ar[i]
    grid = nc.Dataset(tgrid)
    fn = fn_ar[i]
    
    pHr, OmAr, OmCr, pco2r = OmA_3D(grid,carp)

    tdir = '/data/tjarniko/results/hindcast.201905_dayavg_OmA-pH-pCO2/'
    ncname = tdir + fn

    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
    g.createDimension('depths', 40)
    g.createDimension('ydir',898)
    g.createDimension('xdir',398)

    ts = g.createVariable('pH','f4',('depths','ydir','xdir'))
    ts[:] = pHr
    ts2 = g.createVariable('OmA','f4',('depths','ydir','xdir'))
    ts2[:] = OmAr
    ts3 = g.createVariable('OmC','f4',('depths','ydir','xdir'))
    ts3[:] = OmCr
    ts = g.createVariable('pCO2','f4',('ydir','xdir'))
    ts[:] = pco2r[0,:,:]

    f.close()