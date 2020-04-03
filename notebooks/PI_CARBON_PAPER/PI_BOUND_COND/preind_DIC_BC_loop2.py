from __future__ import print_function
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import scipy as sp
import pandas as pd
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
import time

year = 2015
startind = 0
endind = 100
def oned_moxy(tsal, ttemp, tdic, tta, pres_atm, depth_this):

    size_box = np.shape(tdic)
    size_0 = size_box[0]
    size_1= size_box[1]


    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    #tdepth = np.zeros_like(tsra)
    tpressure[:] = pres_atm
    tdepth = np.ravel(depth_this)
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504
    ttera_is = gsw.t_from_CT(tsra,ttera,tzero)

    response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tdepth, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    pHr = pH.reshape(size_0,size_1)
    OmAr = OmegaA.reshape(size_0,size_1)
    pco2r = pco2.reshape(size_0,size_1)
    
    return pHr, OmAr, pco2r

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def find_DIC_corresp_to_pco2(tsal, ttemp, tpco2, tta, pres_atm, depth_this):
    
    steps = 10000
    tsal_r = np.zeros([steps])
    tsal_r[:] = tsal
    ttemp_r = np.zeros([steps])
    ttemp_r[:] = ttemp
    tta_r = np.zeros([steps])
    tta_r[:] = tta * 1e-3
    tpres_r = np.zeros([steps])
    tpres_r[:] = pres_atm
    depth_r = np.zeros([steps])
    depth_r[:] = depth_this
    tzero = np.zeros([steps])

    end_d = 2400
    start_d = 600
    intvl = (end_d - start_d)/steps
    tdic_r = np.arange(start_d,end_d-0.1,intvl) * 1e-3
    
    response_tup = mocsy.mvars(temp=ttemp_r, sal=tsal_r, alk=tta_r, dic=tdic_r, 
                       sil=tzero, phos=tzero, patm=tpres_r, depth=depth_r, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup    
    
    diffmat = pco2 - tpco2
    idx, ans = find_nearest( diffmat,0 )
    
    if ans> 2:
        print('Danger, pco2 found >2 uatm from pco2 given')
#     print(idx)
#     print('difference between real pco2 and pco2 from calc. dic: ',ans)
#     print('DIC found this way:', tdic_r[idx]*1e3)
    fin_dic = tdic_r[idx]*1e3
    
    return fin_dic

def co2_from_year(year, co2_rec):
    '''takes a value for a year, converts year to int,
    and finds appropriate co2 values  from pandas lookup table. 
    if year < 1832, value is for year 1832, if year > 2018, value is for year 2018'''

    whole_year = int(year)

    if whole_year >= 2018:
        whole_year = 2018     
        #print('year > 2018, using value for 2018')

    if whole_year <= 1832:
        whole_year = 1832
        #print('year < 1832, using value for 1832')


    match = (co2_rec['YEAR'] == whole_year) 
    atmco2 = co2_rec['PPMCO2'][match]
    t_co2 = atmco2.values[0]
    return t_co2
    
def preind_dic_ncmaker(startind, endind, year):
#1 open given boundary conditions file and findpco2 and potential density 
    daymon = [31,28,31,30,31,30,31,31,30,31,30,31]
    daymon_LY = [31,29,31,30,31,30,31,31,30,31,30,31]

    year_ar = []
    noday = 365
    if year == 2016:
        t_daymon = daymon_LY
        noday = 366
    else:
        t_daymon = daymon

    for m in range(1,13):
        if m>=10:
            tm = str(m)
        if m<10:
            tm = '0' + str(m)
        print(tm)
        for d in range(1,t_daymon[m-1]+1):
            if d>=10:
                td = str(d)
            if d<10:
                td = '0' + str(d)

            tstr = 'y' + str(year) + 'm' + tm + 'd' + td
            year_ar.append(tstr)
                            
    for ind in range(startind,endind):
        start = time.time()

        print(year_ar[ind])
        test_LO = nc.Dataset('/results/forcing/LiveOcean/boundary_conditions/LiveOcean_v201905_' + year_ar[ind] +'.nc')


        zlevels = (test_LO['deptht'][:])
        sal = test_LO['vosaline'][0,:,0,:]
        temp = test_LO['votemper'][0,:,0,:]
        DIC = test_LO['DIC'][0,:,0,:]
        TA = test_LO['TA'][0,:,0,:]
        depth_this = np.zeros_like(TA)

        #assign depths
        for i in range(0,950):
            depth_this[:,i] = zlevels

        #find potential density and pco2
        potdens = gsw.sigma0(sal,temp)
        pHr, OmAr, pco2r = oned_moxy(sal, temp, DIC, TA, 1, depth_this)
        pco2_resh = pco2r.reshape(40,950)
        print('found pco2 and potential density')
        #find age and co2 intrusion from pot dens
        co2_rec = pd.read_csv('lawdome_maunaloa.csv') 


        #pycnal age from Tetjana method ventilation timescale
        params0 = 0.5406570488955575
        params1 = 11.595484591464743
        params2 = 0.2153367515219747

        pycnal_last_at_surface = 2019 - params0 *np.exp(-params1*(25.15-potdens))+params2
        pycnal_original_co2 = np.zeros_like(pycnal_last_at_surface)
        np.shape(pycnal_original_co2)
        for i in range(0,40):
            for j in range(0,950):
                ty = pycnal_last_at_surface[i,j]
                tco2 = co2_from_year(ty,co2_rec)
                pycnal_original_co2[i,j] = tco2

        #preindustrial pco2 - corrected
        pycnal_intrusion = pycnal_original_co2 - 284
        preind_pco2 = pco2_resh - pycnal_intrusion

        print('found preindustrial pco2')

        #calculate preindustrial DIC from pco2
        preind_dic = np.zeros_like(DIC)
        preind_dic_r = np.ravel(preind_dic)
        pco2r_preind_r = np.ravel(preind_pco2)
        depth_r = np.ravel(depth_this)
        sal_r = np.ravel(sal)
        temp_r = np.ravel(temp)
        DIC_r = np.ravel(DIC)
        TA_r = np.ravel(TA)
        print('calculating preindustrial DIC')
        for i in range(0,len(depth_r)):
            if i%950 == 0:
                print(i)
            t_dic = find_DIC_corresp_to_pco2(sal_r[i], temp_r[i], pco2r_preind_r[i], TA_r[i], 1, depth_r[i])
            preind_dic_r[i] = t_dic

        preind_dic_fin = preind_dic_r.reshape(40,950)

        f = nc.Dataset('./preind_DIC/LO_corrected_' + year_ar[ind] +'_preind_DIC.nc','w', format='NETCDF4') #'w' stands for write
        g = f.createGroup('preindustrial_DIC')
        #g.createDimension('days', len(NO3_mod))
        g.createDimension('xval', 950)
        g.createDimension('depth', 40)
        ts = g.createVariable('preind_dic','f4',('depth','xval'))
        ts[:] = preind_dic_fin

        f.close()
        end = time.time()
        print(' ')
        print('time taken in s')
        print(end - start)