def find_nearest(array, value):
    
    import numpy as np
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

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

def co2_from_year(year):
    import pandas as pd
    '''takes a value for a year, converts year to int,
    and finds appropriate co2 values  from pandas lookup table. 
    if year < 1832, value is for year 1832, if year > 2018, value is for year 2018'''
    co2_rec = pd.read_csv('lawdome_maunaloa.csv') 
    whole_year = int(year)
    
    if whole_year >= 2018:
        whole_year = 2018     
        #('year > 2018, using value for 2018')
        
    if whole_year <= 1832:
        whole_year = 1832
        #('year < 1832, using value for 1832')

    match = (co2_rec['YEAR'] == whole_year) 
    atmco2 = co2_rec['PPMCO2'][match]
    t_co2 = atmco2.values[0]
    return t_co2


def calc_preind_co2_AOU_method(ncname, datestr):
    import numpy as np
    import netCDF4 as nc
    import gsw
    
    """
    - Opens a LiveOcean netcdf file
    - Calculates potential density
    - Uses a derived relationship between potential density (sigma0) and 
    cfc-freon11 estimated watermass age to determine what atmospheric co2
    the watermass saw when it was last at surface.
   - see comments below
    
    Parameters
    ----------
    ncname : name of path + filename of LiveOcean boundary condition file
    for SalishSeaCast, containing temperature, salinity, DIC, TA
    datestr : a string of form y2018m06d01 to write in resulting netcdf name.
    
    Returns
    -------
    writes netcdf file of form ./preind_DIC/LO_intrusion_' + datestr +'_preind_DIC.nc
    with the following variables:
    sigma0, pycnal_last_at_surface, pycnal_witnessed_atm_co2
    insitu_pco2, preind_pco2, preind_dic
    """    
    print(datestr)
    #open dataset & retrieve relevant variables, calculate potential density
    test_LO = nc.Dataset(ncname)
    zlevels = (test_LO['deptht'][:])
    sal = test_LO['vosaline'][0,:,0,:]
    temp = test_LO['votemper'][0,:,0,:]
    sigma0 = gsw.sigma0(sal,temp)
    DIC = test_LO['DIC'][0,:,0,:]
    TA = test_LO['TA'][0,:,0,:]
    O2 = test_LO['OXY'][0,:,0,:]
    depth_this = np.zeros_like(TA)
    #depth_this - array of depths of same shape as DIC
    for i in range(0,950):
        depth_this[:,i] = zlevels
    
    #calculate pycnal's last surfacing, according to exp function
    #found using cfc ages
    params0 = 0.1301889490932413
    params1 = 3.8509914822057825
    params2 = 8.301166081413104
    pycnal_last_at_surface = 2019 - (params0 *np.exp(-params1*(25.15-sigma0))+params2)

    #find last seen atmospheric co2
    pycnal_witnessed_atm_co2 = np.zeros_like(pycnal_last_at_surface)
    for i in range(0,40):
        for j in range(0,950):
            ty = pycnal_last_at_surface[i,j]
            tco2 = co2_from_year(ty)
            pycnal_witnessed_atm_co2[i,j] = tco2
    
#(1) estimate AOU on 26 (assoc with water parcel with DIC_{w,2019,26,jdf})
# = f(O2_{w,2019,26,jdf},S_{w,2019,26,jdf},T_{w,2019,26,jdf}, P_{w,2019,26,jdf})
#(P is there to determine T when last at surface - I'll call it preT next)

    osol = gsw.O2sol(sal,temp,depth_this,-125,50)
    AOU = osol - O2
    print('max AOU: '+str(np.max(AOU)) + ', min AOU: '+ str(np.min(AOU)))
    AOU_stoich = np.copy(AOU)
    AOU_stoich = AOU_stoich * (117/170)
    #AOU_zeroed[AOU<0] = 0         
            
#(2) estimate preformed DIC on 26 when last at surface (say 16 years ago or in 2003): 
# preDIC_{w,2003,26} = DIC_{w,2019,26,jdf} - AOU_{w,2019,26,jdf}
#{AOU may be about 130 umol/kg or so in this e.g.  - we are taking away the AOU because the water hasn't received the organic rain yet - here I am assuming that one
# unit of AOU means one unit of DIC and it will be close but we should check Gruber et al. 2006}
    preformed_DIC = DIC - AOU_stoich

#(3)  estimate preformed PCO2
#prePCO2_{w,jdf,26} = f(preDIC_{w,2003,26},TA_{w,2019,26,jdf},S_{w,2019,26,jdf},preT_{w,2019,26,jdf})
    print('finding preformed_pco2 at surface')
    zeros_here = np.zeros_like(depth_this)
    pHr, OmAr, pco2r = oned_moxy(sal, temp, preformed_DIC, TA, 1, zeros_here)
    preformed_pco2 = pco2r.reshape(40,950)
    print('max preformed_pco2: '+str(np.max(preformed_pco2)) + ', min preformed_pco2: '+ str(np.min(preformed_pco2)))

#(4) estimate disequilibrium PCO2 when last at surface
#diseqPCO2 = prePCO2_{w,jdf,26} - PCO2_a,2003
#{expect diseqPCO2 may be about 0-30 uatm but not 300uatm or more}
    diseqPCO2 = preformed_pco2 - pycnal_witnessed_atm_co2
    print('max diseqPCO2: '+str(np.max(diseqPCO2)) + ', min diseqPCO2: '+ str(np.min(diseqPCO2)))
    pref_pco2_inc_diseqpco2 = diseqPCO2 + 284

    
#(5) estimate preindustrial preformed DIC
#preDIC_PI = f( (PCO2_PI + diseqPCO2), TA_{w,2019,26,jdf},S_{w,2019,26,jdf},preT_{w,2019,26,jdf})

    #pref_pco2_inc_diseqpco2 = diseqPCO2 + preformed_pco2
    print('calculating preindustrial preformed DIC')    
    preind_dic = np.zeros_like(DIC)
    preind_dic_r = np.ravel(preind_dic)
    pref_pco2_inc_diseqpco2_r = np.ravel(pref_pco2_inc_diseqpco2)
    depth_r = np.ravel(depth_this)
    sal_r = np.ravel(sal)
    temp_r = np.ravel(temp)
    TA_r = np.ravel(TA)

    #for i in range(0,len(depth_r)):
    for i in range(0,len(depth_r)):
        if i%950 == 0:
            print(i)
        t_dic = find_DIC_corresp_to_pco2(sal_r[i], temp_r[i], pref_pco2_inc_diseqpco2_r[i], TA_r[i], 1, depth_r[i])
        preind_dic_r[i] = t_dic
    preind_pref_dic = preind_dic_r.reshape(40,950)
    
    deltaDIC = preformed_DIC - preind_pref_dic
    print('max deltaDIC: '+str(np.max(deltaDIC)) + ', min deltaDIC: '+ str(np.min(deltaDIC)))

    final_preind_DIC = DIC - deltaDIC

    f = nc.Dataset('./preind_DIC/LO_AOUmethod_stoicCO_diseq_' + datestr +'_preind_DIC.nc','w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('preindustrial_DIC')
    g.createDimension('xval', 950)
    g.createDimension('depth', 40)
    ts = g.createVariable('sigma0','f4',('depth','xval'))
    ts[:] = sigma0
    ts2 = g.createVariable('pycnal_last_at_surface','f4',('depth','xval'))
    ts2[:] = pycnal_last_at_surface
    ts3 = g.createVariable('pycnal_witnessed_atm_co2','f4',('depth','xval'))
    ts3[:] = pycnal_witnessed_atm_co2
    ts4 = g.createVariable('AOU','f4',('depth','xval'))
    ts4[:] = AOU
    ts5 = g.createVariable('preformed_pco2','f4',('depth','xval'))
    ts5[:] = preformed_pco2
    ts6 = g.createVariable('deltaDIC','f4',('depth','xval'))
    ts6[:] = deltaDIC
    ts7 = g.createVariable('preind_dic','f4',('depth','xval'))
    ts7[:] = final_preind_DIC
    f.close()

def preind_dic_ncmaker(startind, endind, year):
    import time
    import netCDF4 as nc
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
        test_LO = '/results/forcing/LiveOcean/boundary_conditions/LiveOcean_v201905_' + year_ar[ind] +'.nc'
        calc_preind_co2_AOU_method(test_LO, year_ar[ind])