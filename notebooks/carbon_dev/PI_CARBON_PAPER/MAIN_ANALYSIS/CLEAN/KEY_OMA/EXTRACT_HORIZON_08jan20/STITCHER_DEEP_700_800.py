j_st = 700
j_en = 800
shallow = False

#import things
import netCDF4 as nc
import numpy as np
#

##Algorithms

def find_depth_deepalg(dp,prof,water_depth):
    #finds saturation horizon given a profile and corresponding depths
    first_proper_undersat = np.nan
    depth_undersat = np.nan    
    dummy_var = 0
    #print(prof)
    #print('')
    if np.ma.min(prof) >=1e19:
        dummy_var = 0
        depth_undersat = np.nan
        #print('masks all around!')
    elif np.ma.min(prof) >=1:
        depth_undersat = water_depth
        #print('saturated column')
    elif np.ma.max(prof) <1:
        depth_undersat = 0
        #print('undersat to surface')        
    else:
        t_ind = np.where(prof<1)
        t_indar = t_ind[0][0]
        t_indss = np.where(prof>=1)
        t_indsssar = t_indss[0][0]
        if t_indar.size == 0:
            dummy_var = 0
        else:
            if (t_indar.size != 0) & (t_indsssar.size == 0):
                depth_undersat = 0
                first_proper_undersat = 0
                dummy_var = 0
                #print('undersat to surface!')
                max_supsat = np.nan
            else:    
                max_supsat = np.max(t_indsssar)    
                try:
                    first_proper_undersat = np.min(t_indar[t_indar>max_supsat])
                except:
                    dummy_var = 0
                    #print("An exception occurred")
                if first_proper_undersat == 0:
                    depth_undersat = dp[0]
                elif np.isnan(first_proper_undersat):
                    dummy_var = 0
                    #print('saturated watercolumn!')
                else:
                    depth_undersat = (dp[first_proper_undersat]+dp[first_proper_undersat-1])/2
    return depth_undersat

def find_depth_shallowalg(dp,prof,water_depth):
    #finds saturation horizon given a profile and corresponding depths
    first_proper_undersat = np.nan
    depth_undersat = np.nan    
    dummy_var = 0
    #tot masked
    #print(prof)
    #print('')
    if np.ma.min(prof) >=1e19:
        dummy_var = 0
        depth_undersat = np.nan
        #print('masks all around!')
    elif np.ma.min(prof) >=1:
        depth_undersat = water_depth
        #print('saturated column')
    elif np.ma.max(prof) <1:
        depth_undersat = 0
        #print('undersat to surface')        
    else:
        t_ind = np.where(prof<1)
        t_indar = t_ind[0][0]
        t_indss = np.where(prof>=1)
        t_indsssar = t_indss[0][0]
        if t_indar.size == 0:
            dummy_var = 0
        else:
            if (t_indar.size != 0) & (t_indsssar.size == 0):
                depth_undersat = 0
                first_proper_undersat = 0
                dummy_var = 0
                #print('undersat to surface!')
                max_supsat = np.nan
            else:    
                max_supsat = np.max(t_indsssar)    
                try:
                    first_proper_undersat = np.min(t_indar)
                except:
                    dummy_var = 0
                    #print("An exception occurred")
                if first_proper_undersat == 0:
                    depth_undersat = dp[0]
                elif np.isnan(first_proper_undersat):
                    dummy_var = 0
                    #print('saturated watercolumn!')
                else:
                    depth_undersat = (dp[first_proper_undersat]+dp[first_proper_undersat-1])/2
    return depth_undersat




def find_sathor(j_st,j_en,tit,ncname,shallow):
    if shallow == False:
        print('using deep alg')
    if shallow == True:
        print('using shallow alg')
    OmA = nc.Dataset(ncname)
    BR_omA = OmA['model_output']['OmA']


    t_nc = nc.Dataset('/results2/SalishSea/nowcast-green.201812/01jan18/SalishSea_1h_20180101_20180101_grid_T.nc')
    bath =  nc.Dataset('/home/sallen/MEOPAR/grid/bathymetry_201702.nc')
    zlevels = (t_nc['deptht'][:])

    oma_d_BR = np.zeros([365,898,398])
    for day in range(0,365):
        print('ncname is ' + ncname)
        print('day is ' + str(day))
        for j in range(j_st,j_en):
            if j%20 == 0:
                print('j is ' + str(j))
            for i in range(0,398):
                OmA_BR_test = BR_omA[day,:,j,i]
                
                water_depth = bath.variables['Bathymetry'][j,i]
                if shallow == False:
                    oma_dep_BR = find_depth_deepalg(zlevels,OmA_BR_test,water_depth)

                if shallow == True:
                    oma_dep_BR = find_depth_shallowalg(zlevels,OmA_BR_test,water_depth)
                oma_d_BR[day,j,i] = oma_dep_BR

                    
    ncname = tit + str(j_st) + '_' + str(j_en) + '.nc'
    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
#g.createDimension('days', len(3))
    g.createDimension('days', 365)
    g.createDimension('ys', 898)
    g.createDimension('xs', 398)
    ts = g.createVariable('OmAr_HORIZON','f4',('days','ys','xs'))
    ts[:] = oma_d_BR

    f.close()

##CALLS


tit = 'BR_DEEPALG_'
ncname = 'BR_OmA.nc'
find_sathor(j_st,j_en,tit,ncname,shallow)

tit = 'PI_DEEPALG_'
ncname = 'PI_OmA.nc'
find_sathor(j_st,j_en,tit,ncname,shallow)

tit = 'LA_DEEPALG_'
ncname = 'LA_OmA.nc'
find_sathor(j_st,j_en,tit,ncname,shallow)
