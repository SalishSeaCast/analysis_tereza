from numpy import *
from scipy import *
import netCDF4 as nc
import numpy as np
import scipy as sp
import seawater
import datetime as dt
import depth_algs as da

def find_sathor(j_st,j_en,tit,shallow):
    if shallow == False:
        print('using deep alg')
    if shallow == True:
        print('using shallow alg')
    OmA = nc.Dataset('Oma_2015_fixed.nc')
    BR_omA = OmA['model_output']['OmAr_br']
    PI_omA = OmA['model_output']['OmAr_pi']

    t_nc = nc.Dataset('/results2/SalishSea/nowcast-green.201806/01jan18/SalishSea_1h_20180101_20180101_grid_T.nc')
    bath =  nc.Dataset('/home/sallen/MEOPAR/grid/bathymetry_201702.nc')
    zlevels = (t_nc['deptht'][:])
    
    oma_d_PI = np.zeros([365,898,398])
    oma_d_BR = np.zeros([365,898,398])
    for day in range(0,365):
        print(day)
        for j in range(j_st,j_en):
            if j%10 == 0:
                print(j)
            for i in range(0,398):
                OmA_BR_test = BR_omA[day,:,j,i]
                OmA_PI_test = PI_omA[day,:,j,i]
                
                water_depth = bath.variables['Bathymetry'][j,i]
                if shallow == False:
                    oma_dep_BR = da.find_depth_deepalg(zlevels,OmA_BR_test,water_depth)
                    oma_dep_PI = da.find_depth_deepalg(zlevels,OmA_PI_test,water_depth)
                if shallow == True:
                    oma_dep_BR = da.find_depth_shallowalg(zlevels,OmA_BR_test,water_depth)
                    oma_dep_PI = da.find_depth_shallowalg(zlevels,OmA_PI_test,water_depth)
                oma_d_PI[day,j,i] = oma_dep_PI
                oma_d_BR[day,j,i] = oma_dep_BR
                if (i == 250) & (j%10 == 0):
                    print('BR depth: '+str(oma_dep_BR), 'PI depth: '+str(oma_dep_PI)) 
                    
    ncname = tit + str(j_st) + '_' + str(j_en) + '.nc'
    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
#g.createDimension('days', len(3))
    g.createDimension('days', 365)
    g.createDimension('ys', 898)
    g.createDimension('xs', 398)
    ts = g.createVariable('OmArHORIZON_pi','f4',('days','ys','xs'))
    ts[:] = oma_d_PI
    ts2 = g.createVariable('OmArHORIZON_br','f4',('days','ys','xs'))
    ts2[:] = oma_d_BR
    f.close()