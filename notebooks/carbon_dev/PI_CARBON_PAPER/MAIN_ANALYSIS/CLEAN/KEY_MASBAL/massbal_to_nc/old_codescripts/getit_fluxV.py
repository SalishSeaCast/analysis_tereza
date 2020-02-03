import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import scipy as sp
import datetime as dt
""
from salishsea_tools import (
    nc_tools,
    viz_tools,
    geo_tools,
    tidetools
)
import netCDF4 as nc
import cmocean as cm
import glob
import arrow
import gsw

def make_nclen(start,end,ftype, sdir):
    base_ar = []
    sens_ar = []
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
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/MAIN/' + sdir + '/ncs/SKOG_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        #print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        
    return sens_ar
    

start = '2015-01-01'
end = '2015-12-31'
sdir = 'BR_2nd_2015'

ftype = 'dian_V'
varname = 'DIC_VT'

sens_ar = make_nclen(start,end,ftype, sdir)


grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
#get udir, vdir, wdir
vdir = grid['e2t'][0,:,:]
udir = grid['e1t'][0,:,:]
wdir = grid['e3t_0'][0,:,:,:]
tmask = grid['tmask'][0,:,:,:]

#ok so here our tmask is 1 iff we have an ocean, 
print(np.shape(tmask))
print('got this far')
#we set it to 0 at the 20 cell border to not worry aboVT the sides of the domain
#we also set tdat (our variable) that way 
# tmask[:,878:898,:] = 0
# tmask[:,:,0:20] = 0
# tdat[:,878:898,:] = 0
# tdat[:,:,0:20] = 0

#surface area for every cell
surfa = udir*vdir

#broadcast surface area 
surfa_broadcast = np.zeros((40, 898, 398))

for d in range(0,40):
    surfa_broadcast[d,:,:] = surfa[:,:]

#finally we have the cell size for every cell in the system. it's 0 if the cell is not water in the domain we want
cellsize = surfa_broadcast*wdir*tmask
    
    
BR_sums_perday = np.zeros((40,365))
PI_sums_perday = np.zeros((40,365))
LA_sums_perday = np.zeros((40,365))

start = '2015-01-01'
end = '2015-12-31'
ftype = 'dian_V'
varname = 'DIC_VT'
sdir = 'BR_2nd_2015'


BR_ar = make_nclen(start,end,ftype, sdir)

sdir = 'PI_3rd_2015'
LA_ar = make_nclen(start,end,ftype, sdir)

sdir = 'PI_ACBC_2015_3'
PI_ar = make_nclen(start,end,ftype, sdir)

print('done making nclen for dic')

for i in range(0,365):
    
    if (i%5 ==0):
        print(i)
    
    t_test = nc.Dataset(BR_ar[i])
    tdat = t_test[varname][:]
    tdat[:,878:898,:] = 0
    tdat[:,:,0:20] = 0
    tdat_fc = tdat[0,:,:,:]
    tdat_withvol = tdat_fc*cellsize
    tdat_alldomain = np.sum(np.sum(tdat_withvol,axis = 1),axis = 1)

    BR_sums_perday[:,i] =  tdat_alldomain
   
   
import pickle
pickle.dump(BR_sums_perday, open("BR_sums_perday_fV.pkl", 'wb'))
pickle.dump(cellsize, open("cellsize.pkl", 'wb'))

for i in range(0,365):
    
    if (i%15 ==0):
        print(i)
    
    t_test = nc.Dataset(LA_ar[i])
    tdat = t_test[varname][:]
    tdat[:,878:898,:] = 0
    tdat[:,:,0:20] = 0
    tdat_fc = tdat[0,:,:,:]
    tdat_withvol = tdat_fc*cellsize
    tdat_alldomain = np.sum(np.sum(tdat_withvol,axis = 1),axis = 1)

    LA_sums_perday[:,i] =  tdat_alldomain

    
pickle.dump(LA_sums_perday, open("LA_sums_perday_fV.pkl", 'wb'))

for i in range(0,365):
    
    if (i%15 ==0):
        print(i)
    
    t_test = nc.Dataset(PI_ar[i])
    tdat = t_test[varname][:]
    tdat[:,878:898,:] = 0
    tdat[:,:,0:20] = 0
    tdat_fc = tdat[0,:,:,:]
    tdat_withvol = tdat_fc*cellsize
    tdat_alldomain = np.sum(np.sum(tdat_withvol,axis = 1),axis = 1)

    PI_sums_perday[:,i] =  tdat_alldomain

    
pickle.dump(PI_sums_perday, open("PI_sums_perday_fV.pkl", 'wb'))


