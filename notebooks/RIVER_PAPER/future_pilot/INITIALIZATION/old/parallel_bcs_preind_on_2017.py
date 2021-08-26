import netCDF4 as nc
import numpy as np
import pandas as pd
import pickle
import cmocean as cm
import glob
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
import arrow
import gsw
from multiprocessing import Process
import time
import BC_helper_methods as hm

def calc_preind_co2_AOU_method(arrowdate, obs_year, target_year, scen):
    import numpy as np
    import netCDF4 as nc
    import gsw
    
    test_LO = hm.load_nc(arrowdate)
    
    tdate = arrowdate
    yy = tdate.format('YYYY')
    mm = tdate.format('MM')
    dd = tdate.format('DD')
    ymd = f'y{yy}m{mm}d{dd}'
    #open dataset & retrieve relevant variables, calculate potential density

    zlevels = (test_LO['deptht'][:])
    sal = test_LO['vosaline'][0,:,0,:]
    temp = test_LO['votemper'][0,:,0,:]
    sigma0 = gsw.sigma0(sal,temp)
    DIC = test_LO['DIC'][0,:,0,:]
    TA = test_LO['TA'][0,:,0,:]
    O2 = test_LO['OXY'][0,:,0,:]
    depth_this = np.zeros_like(TA)
    zeros = np.zeros_like(TA)
    #depth_this - array of depths of same shape as DIC
    for i in range(0,950):
        depth_this[:,i] = zlevels
 
 ### GET AGE AND WATERMASS WITNESSED CO2   
    #calculate pycnal's last surfacing, according to exp function
    #found using cfc ages
    params0 = 0.1301889490932413
    params1 = 3.8509914822057825
    params2 = 8.301166081413104 #change to 2015 since model year is 2015


    water_age = (params0 *np.exp(-params1*(25.15-sigma0))+params2)
    # year_watermass_at_surface = int(targetyear - age)
    # watermass_witnessed_co2_obs = int(hm.co2_from_year(scen,year_watermass_at_surface))

    #get witnessed co2 both of the present day water parcel and the target year water parcel
    obs_year_ar = np.zeros_like(water_age)
    obs_year_ar[:] = obs_year
    
    target_year_ar = np.zeros_like(water_age)
    target_year_ar[:] = target_year
    
    year_watermass_at_surface = (obs_year_ar - water_age).astype(int)
    watermass_witnessed_co2_obs = hm.co2_from_year(scen,year_watermass_at_surface)
    print(np.min(watermass_witnessed_co2_obs))
    print(np.shape(watermass_witnessed_co2_obs))
#     watermass_witnessed_co2_target = \
#     hm.co2_from_year(scen,year_watermass_at_surface+(target_year_ar-obs_year_ar))
    if target_year < 1905:
        watermass_witnessed_co2_target = 284
    else:
        watermass_witnessed_co2_target = \
    hm.co2_from_year(scen,year_watermass_at_surface+(target_year_ar-obs_year_ar))ttm
###GET AOU
#(1) estimate AOU on 26 (assoc with water parcel with DIC_{w,2019,26,jdf})
# = f(O2_{w,2019,26,jdf},S_{w,2019,26,jdf},T_{w,2019,26,jdf}, P_{w,2019,26,jdf})
#(P is there to determine T when last at surface - I'll call it preT next)
    AOU_stoich = hm.get_AOU_stoich(sal,temp,O2,sigma0,water_age)

 ### GET PREFORMED DIC        
    obs_preformed_dic = DIC - AOU_stoich

#### get preformed pco2 and target year preformed pco2
    pHr, OmAr, pco2r = hm.oned_moxy(sal, temp, obs_preformed_dic, TA, 1, np.zeros_like(sal))
    obsyear_pref_pco2 = pco2r
    diseqPCO2 = obsyear_pref_pco2 - watermass_witnessed_co2_obs
    targetyear_pref_pco2 = watermass_witnessed_co2_target + diseqPCO2

    print('calculating target year preformed DIC')    
    target_preformed_dic = np.zeros_like(DIC)
    target_preformed_dic_r = np.ravel(target_preformed_dic)
    targetyear_pref_pco2_r = np.ravel(targetyear_pref_pco2)
    depth_r = np.ravel(depth_this)
    sal_r = np.ravel(sal)
    temp_r = np.ravel(temp)
    TA_r = np.ravel(TA)
    zeros_r = np.zeros_like(TA_r)
    sigma0_r = np.ravel(sigma0)
 
    start = time.time()
    
    for i in range(0,len(TA_r)):
        if i%(950*5) == 0:
            print(f'level: {i/(950*5)}')
        ### the surface needs better handling    
        if sigma0_r[i] < 25.0:
            target_preformed_dic_r[i] = 9999
        #### the bottom can be actually taken from the cell right above it, no need to do this painful calculation
        if depth_r[i] > 330:
            target_preformed_dic_r[i] = 6666
        else: 
            t_dic = hm.find_DIC_corresp_to_pco2(sal_r[i], temp_r[i], targetyear_pref_pco2_r[i], TA_r[i], 1, 0)
            target_preformed_dic_r[i] = t_dic
    
    
    print('seconds taken at the hard part')
    print(time.time()-start)
    
#     deltaDIC = obs_preformed_dic - target_preformed_dic
            
#     print('max deltaDIC: '+str(np.max(deltaDIC)) + ', min deltaDIC: '+ str(np.min(deltaDIC)))

#     final_target_DIC = DIC - deltaDIC
    
    ## the top and bottom can be dealt with differently

    DIC_r = np.ravel(DIC)
    obs_preformed_dic_r = np.ravel(obs_preformed_dic)
    deltaDIC_r = obs_preformed_dic_r - target_preformed_dic_r
    final_target_DIC_r = DIC_r - deltaDIC_r
    
    for i in range(0,len(TA_r)):
        if sigma0_r[i] < 25.0:
            deltaDIC_r[i] = 9999
            obs_preformed_dic_r[i] = 9999
            target_preformed_dic_r[i] = 9999
            final_target_DIC_r[i] = 9999
            
        if depth_r[i] > 330:    
            deltaDIC_r[i] = 6666
            obs_preformed_dic_r[i] = 6666
            target_preformed_dic_r[i] = 6666
            final_target_DIC_r[i] = 6666

    deltaDIC = deltaDIC_r.reshape(40,950)
    obs_preformed_dic = obs_preformed_dic_r.reshape(40,950)
    target_preformed_dic = target_preformed_dic_r.reshape(40,950)
    final_target_DIC = final_target_DIC_r.reshape(40,950)
    #target_preformed_dic = target_preformed_dic_r.reshape(40,950)
            
#     target_year_ar = np.zeros_like(final_target_DIC)
#     target_year_ar[:] = target_year
    
    f = nc.Dataset(f'./JdF_future_DIC/LO_TY_{target_year}_scen_{scen}_{ymd}_DIC_nosurfnodeep.nc','w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('preindustrial_DIC')
    g.createDimension('xval', 950)
    g.createDimension('depth', 40)
    g.createDimension('single', 1)
    
    ts = g.createVariable('sigma0','f4',('depth','xval'))
    ts[:] = sigma0

    ts2 = g.createVariable('water_age','f4',('depth','xval'))
    ts2[:] = water_age

    ts2a = g.createVariable('target_year','f4',('single'))
    ts2a[:] = target_year  

#     ts3 = g.createVariable('watermass_witnessed_co2_obs','f4',('depth','xval'))
#     ts3[:] = watermass_witnessed_co2_obs

#     ts3a = g.createVariable('watermass_witnessed_co2_target','f4',('depth','xval'))
#     ts3a[:] = watermass_witnessed_co2_target
    
    ts4 = g.createVariable('AOU_stoich','f4',('depth','xval'))
    ts4[:] = AOU_stoich
    
    ts5 = g.createVariable('obsyear_pref_pco2','f4',('depth','xval'))
    ts5[:] = obsyear_pref_pco2
    ts5a = g.createVariable('targetyear_pref_pco2','f4',('depth','xval'))
    ts5a[:] = targetyear_pref_pco2
    
    ts5b = g.createVariable('obsyear_pref_dic','f4',('depth','xval'))
    ts5b[:] = obs_preformed_dic
    ts5c = g.createVariable('targetyear_pref_dic','f4',('depth','xval'))
    ts5c[:] = target_preformed_dic
    
    ts6 = g.createVariable('final_target_DIC','f4',('depth','xval'))
    ts6[:] = final_target_DIC

    f.close()

    
#####pARALELALALALAALA

from multiprocessing import Process

start ='2017-01-01'
end ='2017-12-31'

start_run = arrow.get(start)
end_run = arrow.get(end)

arrow_array = []

for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array.append(r)
    
    
def func1():
  print('func1: starting')
  for i in range(0,60):
        calc_preind_co2_AOU_method(arrow_array[i][0], 2017, 1901, '2_4pt5')

def func2():
  print('func2: starting')
  for i in range(60,120):
        calc_preind_co2_AOU_method(arrow_array[i][0], 2017, 1901, '2_4pt5')
      
def func3():
  print('func3: starting')
  for i in range(120,180):
        calc_preind_co2_AOU_method(arrow_array[i][0], 2017, 1901, '2_4pt5')
        
def func4():
  print('func4: starting')
  for i in range(180,240):
        calc_preind_co2_AOU_method(arrow_array[i][0], 2017, 1901, '2_4pt5')
        
def func5():
  print('func5: starting')
  for i in range(240,300):
        calc_preind_co2_AOU_method(arrow_array[i][0], 2017, 1901, '2_4pt5')
        
def func6():
  print('func36: starting')
  for i in range(300,366):
        calc_preind_co2_AOU_method(arrow_array[i][0], 2017, 1901, '2_4pt5')
        
if __name__ == '__main__':
  p1 = Process(target=func1)
  p1.start()
  p2 = Process(target=func2)
  p2.start()
  p3 = Process(target=func3)
  p3.start()
  p4 = Process(target=func4)
  p4.start()
  p5 = Process(target=func5)
  p5.start()
  p6 = Process(target=func6)
  p6.start()                
    
  p1.join()
  p2.join()
  p3.join()
  p4.join()
  p5.join()
  p6.join()