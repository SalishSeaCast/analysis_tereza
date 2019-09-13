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


def calc_preind_co2_intrusion_method(ncname, datestr):
    import numpy as np
    import netCDF4 as nc
    import gsw
    
    """
    - Opens a LiveOcean netcdf file
    - Calculates potential density
    - Uses a derived relationship between potential density (sigma0) and 
    cfc-freon11 estimated watermass age to determine what atmospheric co2
    the watermass saw when it was last at surface. 
    - Then calculates pco2 intrusion = (witnessed atm co2 - preindustrial atm co2)
    - Then calculates in-situ present-day pco2 using mocsy
    - Then subtracts intrusion from present-day pco2 to get preindustrial pco2
    - Then iteratively finds preindustrial DIC corresponding to preindustrial pco2
    - Saves these variables in './preind_DIC/LO_intrusion_' + datestr +'_preind_DIC.nc' 
    
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
            
    #calculate pco2 intrusion - how much 'extra' co2 did this pycnal
    #see, compared to preindustrial levels?
    pycnal_intrusion = pycnal_witnessed_atm_co2 - 284
    
    print('calculate parcel in-situ pco2')
    #calculate parcel's in-situ pco2
    pHr, OmAr, pco2r = oned_moxy(sal, temp, DIC, TA, 1, depth_this)
    insitu_pco2 = pco2r.reshape(40,950)
    
    #calculate preindustrial pco2 = insitu_pco2 - intrusion
    preind_pco2 = insitu_pco2 - pycnal_intrusion
    
    #calculate preindustrial DIC iteratively
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
    preind_dic = preind_dic_r.reshape(40,950)
    
    f = nc.Dataset('./preind_DIC/LO_intrusion_' + datestr +'_preind_DIC.nc','w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('preindustrial_DIC')
    g.createDimension('xval', 950)
    g.createDimension('depth', 40)
    ts = g.createVariable('sigma0','f4',('depth','xval'))
    ts[:] = sigma0
    ts2 = g.createVariable('pycnal_last_at_surface','f4',('depth','xval'))
    ts2[:] = pycnal_last_at_surface
    ts3 = g.createVariable('pycnal_witnessed_atm_co2','f4',('depth','xval'))
    ts3[:] = pycnal_witnessed_atm_co2
    ts4 = g.createVariable('insitu_presday_pco2','f4',('depth','xval'))
    ts4[:] = insitu_pco2
    ts5 = g.createVariable('preind_pco2','f4',('depth','xval'))
    ts5[:] = preind_pco2
    ts6 = g.createVariable('preind_dic','f4',('depth','xval'))
    ts6[:] = preind_dic
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
        calc_preind_co2_intrusion_method(test_LO, year_ar[ind])