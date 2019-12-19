#extract_avg_monthly_oct.py
#oct
start = '2015-10-01'
end = '2015-10-31'
nd = 31

from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import pickle

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import scipy as sp
import warnings
warnings.filterwarnings('ignore')

import seawater
import datetime as dt
""
from salishsea_tools import (
    nc_tools,
    viz_tools,
    geo_tools,
    tidetools
)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use('seaborn-whitegrid')
import netCDF4 as nc

import cmocean as cm
import glob
import sys
sys.path.append('/data/tjarniko/mocsy')
sys.path.append('/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/CCCmaDEV/CCCma_src')
import mocsy
import CCCma
import CCCma_stations as cs
from matplotlib import reload
import arrow
import gsw

import xarray as xr
from datetime import datetime
from salishsea_tools import grid_tools, viz_tools

from math import log10, floor

def mocsy_3d_getOmA(tsal,ttemp,tdic,tta):
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
                       sil=tzero, phos=tzero, patm=tpressure, depth=tzero, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup


    OmegaAR = OmegaA.reshape(40,898,398)
    return OmegaAR

def make_nclist(start,end, var, tdir, verbose):
    
    rootdir = '/data/tjarniko/results/BASERUN_EXP/MAIN/'
    nc_array = []
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    
    
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    dayslen = len(arrow_array)
    
    for i in range(0,dayslen):
        tdate = arrow_array[i][0]
        ymd = tdate.format('YYYYMMDD')
        if verbose == True:
            print(ymd)

        ##SKOG_1d_20151127_20151231_ptrc_T_20151231-20151231.nc    
        t_nc = rootdir + tdir + '/ncs/' \
        '/SKOG_1d_*' + '_'+var+'_' +str(ymd) + '-' + str(ymd) + '.nc'
        t_ncname = glob.glob(t_nc)
        if verbose == True:
            print('file is: ' + t_ncname[0])
            print('')
        nc_array.append(t_ncname[0])
        
    return nc_array


############oct#################

##
tdir = 'PI_ACBC_2015_3'
verbose = False
#oct
var ='carp_T'
nc_C_PI_oct = make_nclist(start,end, var, tdir, verbose)
var = 'grid_T'
nc_TS_PI_oct = make_nclist(start,end, var, tdir, verbose)
#
tdir = 'BR_2nd_2015'
verbose = False
#oct
var ='carp_T'
nc_C_BR_oct = make_nclist(start,end, var, tdir, verbose)
var = 'grid_T'
nc_TS_BR_oct = make_nclist(start,end, var, tdir, verbose)

PI_DIC_oct = np.zeros([nd,40,898,398])
PI_TA_oct = np.zeros([nd,40,898,398])
PI_T_oct = np.zeros([nd,40,898,398])
PI_S_oct = np.zeros([nd,40,898,398])
#
BR_DIC_oct = np.zeros([nd,40,898,398])
BR_TA_oct = np.zeros([nd,40,898,398])
BR_T_oct = np.zeros([nd,40,898,398])
BR_S_oct = np.zeros([nd,40,898,398])
#
for i in range(0,len(nc_TS_PI_oct)):
    if (i%5 == 0):
        print(i)
    grid_PI = nc.Dataset(nc_TS_PI_oct[i])
    carp_PI = nc.Dataset(nc_C_PI_oct[i])
    grid_BR = nc.Dataset(nc_TS_BR_oct[i])
    carp_BR = nc.Dataset(nc_C_BR_oct[i])
    #carbon
    PI_DIC_oct[i,:,:,:]=carp_PI['dissolved_inorganic_carbon'][0,:,:,:]
    PI_TA_oct[i,:,:,:]=carp_PI['total_alkalinity'][0,:,:,:]
    BR_DIC_oct[i,:,:,:]=carp_BR['dissolved_inorganic_carbon'][0,:,:,:]
    BR_TA_oct[i,:,:,:]=carp_BR['total_alkalinity'][0,:,:,:]
    #TS
    PI_T_oct[i,:,:,:]=grid_PI['votemper'][0,:,:,:]
    PI_S_oct[i,:,:,:]=grid_PI['vosaline'][0,:,:,:]
    BR_T_oct[i,:,:,:]=grid_BR['votemper'][0,:,:,:]
    BR_S_oct[i,:,:,:]=grid_BR['vosaline'][0,:,:,:]

PI_DIC_oct_mean = np.nanmean(PI_DIC_oct,axis=0)
PI_TA_oct_mean = np.nanmean(PI_TA_oct,axis=0)
PI_T_oct_mean = np.nanmean(PI_T_oct,axis=0)
PI_S_oct_mean = np.nanmean(PI_S_oct,axis=0)

BR_DIC_oct_mean = np.nanmean(BR_DIC_oct,axis=0)
BR_TA_oct_mean = np.nanmean(BR_TA_oct,axis=0)
BR_T_oct_mean = np.nanmean(BR_T_oct,axis=0)
BR_S_oct_mean = np.nanmean(BR_S_oct,axis=0)

pickle.dump(PI_DIC_oct_mean, open("PI_DIC_oct_mean.pkl", 'wb'))
pickle.dump(PI_TA_oct_mean, open("PI_TA_oct_mean.pkl", 'wb'))
pickle.dump(PI_T_oct_mean, open("PI_T_oct_mean.pkl", 'wb'))
pickle.dump(PI_S_oct_mean, open("PI_S_oct_mean.pkl", 'wb'))

pickle.dump(BR_DIC_oct_mean, open("BR_DIC_oct_mean.pkl", 'wb'))
pickle.dump(BR_TA_oct_mean, open("BR_TA_oct_mean.pkl", 'wb'))
pickle.dump(BR_T_oct_mean, open("BR_T_oct_mean.pkl", 'wb'))
pickle.dump(BR_S_oct_mean, open("BR_S_oct_mean.pkl", 'wb'))

PI_OmA_oct_mean = mocsy_3d_getOmA(PI_S_oct_mean,PI_T_oct_mean,\
                                  PI_DIC_oct_mean,PI_TA_oct_mean)

BR_OmA_oct_mean = mocsy_3d_getOmA(BR_S_oct_mean,BR_T_oct_mean,\
                                  BR_DIC_oct_mean,BR_TA_oct_mean)

pickle.dump(PI_OmA_oct_mean, open("PI_OmA_oct_mean.pkl", 'wb'))
pickle.dump(BR_OmA_oct_mean, open("BR_OmA_oct_mean.pkl", 'wb'))

f = open('PI_OmA_oct_mean.pkl', 'rb')
PI_OmA_oct_mean = pickle.load(f)   
f = open('BR_OmA_oct_mean.pkl', 'rb')
BR_OmA_oct_mean = pickle.load(f) 

mesh = xr.open_dataset('/home/sallen/MEOPAR/grid/mesh_mask201702.nc')
mbathys = np.array(mesh.mbathy[0])
mbathys_cor = mbathys -1 

PI_OmA_oct_benthos = np.zeros([898,398])
BR_OmA_oct_benthos = np.zeros([898,398])

mbathys_cor = mbathys-1
for j in range(0,898):
    if j%50 == 0:
        print(j)
    for i in range(0,398):
        benth = int(mbathys_cor[j,i])
        if benth < 0:
            PI_OmA_oct_benthos[j,i] = 1e20
            BR_OmA_oct_benthos[j,i] = 1e20
   
        else:   
            PI_OmA_oct_benthos[j,i] = PI_OmA_oct_mean[benth,j,i]
            BR_OmA_oct_benthos[j,i] = BR_OmA_oct_mean[benth,j,i]

pickle.dump(PI_OmA_oct_benthos, open("PI_OmA_oct_benthos_mbath.pkl", 'wb'))
pickle.dump(BR_OmA_oct_benthos, open("BR_OmA_oct_benthos_mbath.pkl", 'wb'))
