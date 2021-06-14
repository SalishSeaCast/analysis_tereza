def find_nearest(array, value):
    ## find index of closest value in an array to a given value
    ## used in find_DIC_corresp_to_pco2

    import numpy as np
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def find_DIC_corresp_to_pco2(tsal, ttemp, tpco2, tta, pres_atm, depth_this):
    
    ### given a pco2, salinity, temp, etc get DIC
    import numpy as np
    import mocsy
    import gsw
    
    steps = 10000
    tsal_r = np.zeros([steps])
    tsal_r[:] = tsal
    #convert to psu
    tsal_r_psu = tsal_r*35/35.16504
    
    ttemp_r = np.zeros([steps])
    ttemp_r[:] = ttemp
    #convert temperature to potential temperature
    ttemp_r_pot = gsw.pt_from_CT(tsal_r,ttemp_r)
    tta_r = np.zeros([steps])
    tta_r[:] = tta * 1e-3
    tpres_r = np.zeros([steps])
    tpres_r[:] = pres_atm
    depth_r = np.zeros([steps])
    depth_r[:] = depth_this
    tzero = np.zeros([steps])
    tlat = np.zeros([steps])
    tlat[:]=50
    
    end_d = 2400
    start_d = 600
    intvl = (end_d - start_d)/steps
    tdic_r = np.arange(start_d,end_d-0.1,intvl) * 1e-3
    #change to take potential temperature
    response_tup = mocsy.mvars(temp=ttemp_r_pot, sal=tsal_r_psu, alk=tta_r, dic=tdic_r, 
                       sil=tzero, phos=tzero, patm=tpres_r, depth=depth_r, lat=tzero, 
                        optcon='mol/m3', optt='Tpot', optp='m',
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

    ### GET pCO2 (and Omega, etc) GIVEN DIC, TA 
    import sys
    sys.path.append('/data/tjarniko/mocsy')
    import mocsy
    import numpy as np
    import gsw
    
    size_box = np.shape(tdic)
    size_0 = size_box[0]
    size_1= size_box[1]

    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    #convert cons. temperature to potential temperature
    ttera_pot = gsw.pt_from_CT(tsra,ttera)
    
    
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    #tdepth = np.zeros_like(tsra)
    tpressure[:] = pres_atm
    tdepth = np.ravel(depth_this)
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504
    #ttera_is = gsw.t_from_CT(tsra,ttera,tzero)

    response_tup = mocsy.mvars(temp=ttera_pot, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tdepth, lat=tzero, 
                        optcon='mol/m3', optt='Tpot', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    pHr = pH.reshape(size_0,size_1)
    OmAr = OmegaA.reshape(size_0,size_1)
    pco2r = pco2.reshape(size_0,size_1)
    
    return pHr, OmAr, pco2r

def co2_from_year_pointmeas(scen,tyear):
    import pandas as pd
    import numpy as np

    #get witnessed co2, currently only takes 2-4.5 and 5-8.5

    df = pd.read_csv('./meinshausen/meinshausen_scenario_conc.csv')
    data_top = list(df.columns) 
#     print('available data columns')
#     print(data_top)

    year = np.array(df['year'][:])
    # SSP1_1pt9_NH = np.array(df['SSP1-1.9_NH'][:])
    # SSP1_2pt6_NH = np.array(df['SSP1-2.6_NH'][:])
    # SSP2_4pt5_NH = np.array(df['SSP2-4.5'][:])
    SSP2_4pt5_NH = np.array(df['SSP2-4.5'][:])
    SSP5_8pt5_NH = np.array(df['SSP5-8.5'][:])

    q = np.shape(tyear)
    tco2_ar = np.zeros_like(tyear)
    tyear_r = np.ravel(tyear)
    tco2_ar_r = np.ravel(tco2_ar)
    for i in range(0,len(tco2_ar_r)):
        # building for 2_4pt5
        if (scen == '2_4pt5'): 

            tco2_ar_r[i] = SSP2_4pt5_NH[year == (int(tyear_r[i]))]

        if scen == '5_8pt5':

            tco2_ar_r[i] = SSP5_8pt5_NH[year == (int(tyear_r[i]))]
    
    #tco2_ar = tco2_ar_r.reshape(q[0],q[1])

    return tco2_ar


def co2_from_year(scen,tyear):
    import pandas as pd
    import numpy as np

    #get witnessed co2, currently only takes 2-4.5 and 5-8.5

    df = pd.read_csv('./meinshausen/meinshausen_scenario_conc.csv')
    data_top = list(df.columns) 
#     print('available data columns')
#     print(data_top)

    year = np.array(df['year'][:])
    # SSP1_1pt9_NH = np.array(df['SSP1-1.9_NH'][:])
    # SSP1_2pt6_NH = np.array(df['SSP1-2.6_NH'][:])
    # SSP2_4pt5_NH = np.array(df['SSP2-4.5'][:])
    SSP2_4pt5_NH = np.array(df['SSP2-4.5'][:])
    SSP5_8pt5_NH = np.array(df['SSP5-8.5'][:])

    q = np.shape(tyear)
    tco2_ar = np.zeros_like(tyear)
    tyear_r = np.ravel(tyear)
    tco2_ar_r = np.ravel(tco2_ar)
    for i in range(0,len(tco2_ar_r)):
        # building for 2_4pt5
        if (scen == '2_4pt5'): 

            tco2_ar_r[i] = SSP2_4pt5_NH[year == (int(tyear_r[i]))]

        if scen == '5_8pt5':

            tco2_ar_r[i] = SSP5_8pt5_NH[year == (int(tyear_r[i]))]
    
    tco2_ar = tco2_ar_r.reshape(q[0],q[1])

    return tco2_ar

def load_nc(arrowdate):
    import arrow 
    import netCDF4 as nc

    tdate = arrowdate
    print(tdate)
    
    #take a boundary condition string from liveocean, open it
    yy = tdate.format('YYYY')
    mm = tdate.format('MM')
    dd = tdate.format('DD')
    ymd = f'y{yy}m{mm}d{dd}'
    tstr = f'/results/forcing/LiveOcean/boundary_conditions/LiveOcean_v201905_{ymd}.nc'
    
    w = nc.Dataset(tstr)

    return w

def get_AOU_stoich(sal,temp,O2,sigma0,water_age):
    
    import gsw
    import numpy as np
    zeros = np.zeros_like(sal)
    osol = gsw.O2sol(sal,temp,zeros,-125,50)
    #convert osol to umol/L
    osol_umolL = osol*(1000/(1000+sigma0))
    AOU = osol_umolL - O2
    #print('max AOU: '+str(np.max(AOU)) + ', min AOU: '+ str(np. min(AOU)))
    AOU_stoich = np.copy(AOU)
    AOU_stoich = AOU_stoich * (117/170)
    #AOU_zeroed[AOU<0] = 0 
    ## AOU for young waters?!  
    ##AOU_stoich[water_age<10] = 0

    return AOU_stoich