from __future__ import print_function
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
import pickle
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
#from matplotlib import reload
import arrow
import gsw
import time

def co2_from_year(year):
    import pandas as pd
    '''takes a value for a year, converts year to int,
    and finds appropriate co2 values  from pandas lookup table. 
    if year < 1832, value is for year 1832, if year > 2018, value is for year 2018'''
    co2_rec = pd.read_csv('co2_historical_ssp585.csv') 
    whole_year = int(year)
        
    if whole_year <= 1832:
        whole_year = 1832
        #('year < 1832, using value for 1832')

    match = (co2_rec['YEAR'] == whole_year) 
    atmco2 = co2_rec['PPMCO2'][match]
    t_co2 = atmco2.values[0]
    return t_co2

df = pd.read_csv('coastal_glodap_with_age_3.csv')
#df = pd.read_csv('coastal_glodap_with_age.csv')
data_top = list(df.columns) 
print(data_top)

#need - bottom depth, pressure, dic, ta, salinity, temp, lat, lon
tALK = np.array(df['talk'][:])
tALK_orig = np.copy(tALK)
tDIC = np.array(df['tco2'][:])
tDIC_orig = np.copy(tDIC)
tSAL = np.array(df['salinity'][:])
tTEMP = np.array(df['temperature'][:])
tPRES = np.array(df['pressure'][:])
tLAT = np.array(df['latitude'][:])
tLON = np.array(df['longitude'][:])
tBOTdepth = np.array(df['bottomdepth'][:])
tYEAR = np.array(df['year'])
tAOU = np.array(df['aou'])
tAOU_orig = np.array(df['aou'])
#convert from umol/kg to mmol/m3
import seawater
#help(seawater.dens)
dens = seawater.dens(tSAL,tTEMP,tPRES)
tDIC=tDIC*dens/1000
tALK=tALK*dens/1000
tAOU = tAOU*dens/1000
tAGE = np.array(df['age'])
tALK_DIC = tALK-tDIC
tALK_DIC2 = tALK-(tDIC+50)

print(np.shape(tDIC))

#dic, ta, actual reasonable numbers
filt_ALK = ((tALK>-999) & (~np.isnan(tALK)))
filt_DIC = ((tDIC>-999) & (~np.isnan(tDIC)))
filt_SAL = (tSAL >-999) & (~np.isnan(tSAL))
filt_TEMP = (tTEMP >-999) & (~np.isnan(tTEMP))
filt_PRES = (tPRES >-999) & (~np.isnan(tPRES))
filt_EST = (filt_SAL) & (tSAL >= 20 ) & (tPRES < 201)
filt_AOU = (tAOU>-999) & (~np.isnan(tAOU))
#bottom depth relatively shallow <
filt_DEPTH = (tBOTdepth < 1001) & filt_ALK & filt_DIC & filt_SAL
#year is modern
filt_DEPTH2 = (tBOTdepth < 501) & filt_ALK & filt_DIC & filt_SAL
filt_YEAR = (tYEAR > 2000)
filt_AGE = (tAGE > -9999)

filt_comp = filt_DEPTH2& filt_EST & filt_PRES & filt_SAL & \
filt_TEMP & filt_DIC & filt_ALK & filt_AOU & filt_YEAR & filt_AGE




filt_ALK_DIC = (np.abs(tALK_DIC) < 50) & (filt_ALK) & (filt_DIC) & (filt_SAL)
filt_ALK_DIC2 = (np.abs(tALK_DIC2) < 50) & (filt_ALK) & (filt_DIC) & (filt_SAL)

#glodap v2 is 1275558 
print('total datapoints here')
print(np.shape(tDIC))
print('total datapoints that have an age associated')
print(np.shape(np.where(filt_AGE)))

print('total datapoints with our filt-comp')
print(np.shape(np.where(filt_comp)))
print('total datapoints with our filt-comp with abs(TA-DIC)<50')
print(np.shape(np.where(filt_comp&filt_ALK_DIC)))
print('total datapoints with our filt-comp with abs(TA-(DIC+50))<50')
print(np.shape(np.where(filt_comp&filt_ALK_DIC2)))

tALK_coastal = tALK[filt_comp]
tDIC_coastal = tDIC[filt_comp]
tSAL_coastal=tSAL[filt_comp]
tTEMP_coastal=tTEMP[filt_comp]
tPRES_coastal=tPRES[filt_comp]
tLAT_coastal=tLAT[filt_comp]
tLON_coastal=tLON[filt_comp]
tYEAR_coastal=tYEAR[filt_comp]
tAOU_coastal=tAOU[filt_comp]
tAGE_coastal=tAGE[filt_comp]
tYEAR_coastal = tYEAR[filt_comp]
tYEAR_coastal = tYEAR[filt_comp]


def find_DIC_corresp_to_pco2(tsal, ttemp, tpco2, tta, pres_atm, depth_this):
    
    import numpy as np
    import mocsy
    
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

def oned_moxy(tsal, ttemp, tdic, tta, pres_atm, depth_this):
    import sys
    sys.path.append('/data/tjarniko/mocsy')
    import mocsy
    import numpy as np
    import gsw
    
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
    
    #no need to convert these, they're already converted
    #tsra_psu = tsra*35/35.16504
    #ttera_is = gsw.t_from_CT(tsra,ttera,tzero)

    response_tup = mocsy.mvars(temp=ttera, sal=tsra, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tdepth, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup
    
    return pH, OmegaA, pco2

def find_nearest(array, value):
    
    import numpy as np
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def future_co2_fc(obs_year,target_year,water_age,AOU,DIC,TA,S,T,p, verbose = True):
    """
    Usage: adaptation of c* gruber method 
    Inputs:
        current_year 
        age_of_waterage of the water parcel - how old is the water
        DIC  : umol/L
        TA   : umol/L
        AOU  : umol/L
        S    :  g/kg
        T    :  cons. temp
        p    :  meters/dbar
        
    """

    year_watermass_at_surface = obs_year - water_age
    watermass_witnessed_co2_obs = co2_from_year(year_watermass_at_surface)
    watermass_witnessed_co2_target = co2_from_year(year_watermass_at_surface+(target_year-obs_year))
    #convert AOU to AOU_stoich, ratio of 
    AOU_stoich = AOU * (117/170)
    
    #step A
    preformed_DIC = DIC - AOU_stoich
        
    if verbose == True:

        print('DIC of observation: '+str(DIC) + ' umol/kg')
        print('AOU: '+str(AOU) + ' umol/kg')
        print('stoich AOU (mult. by C:O ratio 117/170): '+str(AOU_stoich) + ' umol/kg')    
        print('preformed_DIC (DIC when last surface, so DIC less stoich AOU) '+str(preformed_DIC))
        #print('year to calculate co2 for: '+str(current_year))
        print('watermass age: ' + str(water_age))
        print('year observation was made: '+str(obs_year))
        print('year we are trying to calculate DIC for: '+str(target_year))
    #get preformed pCO2
    #step B - calculate pref_pco2_when_last_at_surface
    pHr, OmAr, pco2r = oned_moxy(S, T, preformed_DIC, TA, 1, 0)
    obsyear_pref_pco2 = pco2r
    diseqPCO2 = obsyear_pref_pco2 - watermass_witnessed_co2_obs
    targetyear_pref_pco2 = watermass_witnessed_co2_target + diseqPCO2
    
    preind_pref_pco2 = 284 + diseqPCO2
    

    if verbose == True:
        print('watermass last at surface: '+str(year_watermass_at_surface))
        print('watermass witnessed atm co2: ' + str(watermass_witnessed_co2_obs))
        print('preformed_pco2): '+str(pco2r))
        print('diseq co2: '+str(diseqPCO2))
        print('preindustrial preformed pco2 (284+diseq): '+str(preind_pref_pco2))

    #step C: get preindustrial preformed DIC and preindustrial dic
    preind_pref_DIC = find_DIC_corresp_to_pco2(S, T, preind_pref_pco2, TA, 1, 0)
    obsyear_pref_DIC = find_DIC_corresp_to_pco2(S, T, obsyear_pref_pco2, TA, 1, 0)
    targetyear_pref_DIC = find_DIC_corresp_to_pco2(S, T, targetyear_pref_pco2, TA, 1, 0)
    
    preind_DIC = preind_pref_DIC + AOU_stoich
    obsyear_DIC = obsyear_pref_DIC + AOU_stoich
    targetyear_DIC = targetyear_pref_DIC + AOU_stoich
    if verbose == True:
        print('preindustrial preformed DIC: '+str(preind_pref_DIC))
        print('preindustrial DIC: '+str(preind_DIC))
        print('obsyear preformed DIC: '+str(obsyear_pref_DIC))
        print('obsyear DIC: '+str(obsyear_DIC))    
        print('targetyear preformed DIC: '+str(targetyear_pref_DIC))
        print('targetyear DIC: '+str(targetyear_DIC))    

    return targetyear_DIC

def future_co2_fc(obs_year,target_year,water_age,AOU,DIC,TA,S,T,p, verbose = True):
    """
    Usage: adaptation of c* gruber method 
    Inputs:
        current_year 
        age_of_waterage of the water parcel - how old is the water
        DIC  : umol/L
        TA   : umol/L
        AOU  : umol/L
        S    :  g/kg
        T    :  cons. temp
        p    :  meters/dbar
        
    """

    year_watermass_at_surface = obs_year - water_age
    watermass_witnessed_co2_obs = co2_from_year(year_watermass_at_surface)
    watermass_witnessed_co2_target = co2_from_year(year_watermass_at_surface+(target_year-obs_year))
    #convert AOU to AOU_stoich, ratio of 
    AOU_stoich = AOU * (117/170)
    
    #step A
    preformed_DIC = DIC - AOU_stoich
        
    if verbose == True:

        print('DIC of observation: '+str(DIC) + ' umol/kg')
        print('AOU: '+str(AOU) + ' umol/kg')
        print('stoich AOU (mult. by C:O ratio 117/170): '+str(AOU_stoich) + ' umol/kg')    
        print('preformed_DIC (DIC when last surface, so DIC less stoich AOU) '+str(preformed_DIC))
        print('year to calculate co2 for: '+str(current_year))
        print('watermass age: ' + str(water_age))
        print('year observation was made: '+str(obs_year))
        print('year we are trying to calculate DIC for: '+str(target_year))
    #get preformed pCO2
    #step B - calculate pref_pco2_when_last_at_surface
    pHr, OmAr, pco2r = oned_moxy(S, T, preformed_DIC, TA, 1, 0)
    obsyear_pref_pco2 = pco2r
    diseqPCO2 = obsyear_pref_pco2 - watermass_witnessed_co2_obs
    targetyear_pref_pco2 = watermass_witnessed_co2_target + diseqPCO2
    
    preind_pref_pco2 = 284 + diseqPCO2
    

    if verbose == True:
        print('watermass last at surface: '+str(year_watermass_at_surface))
        print('watermass witnessed atm co2: ' + str(watermass_witnessed_co2_obs))
        print('preformed_pco2): '+str(pco2r))
        print('diseq co2: '+str(diseqPCO2))
        print('preindustrial preformed pco2 (284+diseq): '+str(preind_pref_pco2))

    #step C: get preindustrial preformed DIC and preindustrial dic
    preind_pref_DIC = find_DIC_corresp_to_pco2(S, T, preind_pref_pco2, TA, 1, 0)
    obsyear_pref_DIC = find_DIC_corresp_to_pco2(S, T, obsyear_pref_pco2, TA, 1, 0)
    targetyear_pref_DIC = find_DIC_corresp_to_pco2(S, T, targetyear_pref_pco2, TA, 1, 0)
    
    preind_DIC = preind_pref_DIC + AOU_stoich
    obsyear_DIC = obsyear_pref_DIC + AOU_stoich
    targetyear_DIC = targetyear_pref_DIC + AOU_stoich
    if verbose == True:
        print('preindustrial preformed DIC: '+str(preind_pref_DIC))
        print('preindustrial DIC: '+str(preind_DIC))
        print('obsyear preformed DIC: '+str(obsyear_pref_DIC))
        print('obsyear DIC: '+str(obsyear_DIC))    
        print('targetyear preformed DIC: '+str(targetyear_pref_DIC))
        print('targetyear DIC: '+str(targetyear_DIC))    

    return targetyear_DIC

years_record = np.arange(1800,2101,1)

lookup_table = np.zeros([np.size(tDIC_coastal),np.size(years_record)])
print(np.shape(lookup_table))

for i in range(0,len(tALK_coastal)):

    print('i is' + str(i))
    
    for j in range(0,len(years_record)):
        if (j%100 == 0):
            print(j)
        obsyear = tYEAR_coastal[i]
        targetyear = years_record[j]
        water_age = tAGE_coastal[i]
        AOU = tAOU_coastal[i]
        DIC = tDIC_coastal[i]
        TA = tALK_coastal[i]
        S = tSAL_coastal[i]
        T = tTEMP_coastal[i]
        p = tPRES_coastal[i]
        
        
        tyd = future_co2_fc(obsyear,targetyear,water_age,AOU,DIC,TA,S,T,p,verbose = False)
        lookup_table[i,j] = tyd

pickle.dump(lookup_table, open("lookup_table_coastaldat3.pkl", 'wb'))