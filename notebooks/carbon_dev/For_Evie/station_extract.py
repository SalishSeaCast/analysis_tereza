import datetime as dt
import subprocess
import numpy as np
import os
import placesEvie as places
import glob
import time
import sys
import arrow
import netCDF4 as nc

start = '2015-01-01'
end = '2015-12-31'

start_run = arrow.get(start)
end_run = arrow.get(end)
arrow_array = []
date_array = []
for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array.append(r)

#print('days: '+str(len(arrow_array)))
dayslen = len(arrow_array)

for i in range(0,dayslen):
    tdate = arrow_array[i][0]
    ddmmmyy = tdate.format('YYYYMMDD').lower()
    date_array.append(ddmmmyy)

plcnames = ['Cortes/Marina','Deep Bay','Fanny Bay','Lasqueti Island','Lund/Desolation Sound',\
           'Main SoG','Maple Bay','Mouth of Okeover','Nanoose Bay','Northern Baynes',\
           'Salt Spring','Southern Baynes']

plcnamesf = ['CortesMarina','DeepBay','FannyBay','LasquetiIsland','LundDesolation Sound',\
           'MainSoG','MapleBay','MouthofOkeover','NanooseBay','NorthernBaynes',\
           'SaltSpring','SouthernBaynes']


for p in range(0,len(plcnames)):
    
    tplc = plcnames[p]
    tplcf = plcnamesf[p]
    print(tplc)
    
    ji = places.PLACES[tplc]['NEMO grid ji']
    tj = ji[0]
    ti = ji[1]

    TA_stor = np.zeros([365,40])
    DIC_stor = np.zeros([365,40])
    SAL_stor = np.zeros([365,40])
    TEMP_stor = np.zeros([365,40])

    for i in range(0,len(date_array)):
        ymd = date_array[i]
        if i % 50 == 0:
            print(i)
        sdir = '/data/tjarniko/results/BR_1st_2015/ncs/'
        carp = sdir + '/SKOG_1d_*'+ 'carp_T_' + ymd + '-' + ymd + '.nc'
        grid = sdir + '/SKOG_1d_*'+ 'grid_T_' + ymd + '-' + ymd + '.nc'
        tcarp = glob.glob(carp)
        tgrid = glob.glob(grid)
        carpnc = nc.Dataset(tcarp[0])
        gridnc = nc.Dataset(tgrid[0])
        TA = carpnc.variables['total_alkalinity']
        DIC = carpnc.variables['dissolved_inorganic_carbon']
        TEMP = gridnc.variables['votemper']
        SAL = gridnc.variables['vosaline']

        loc_ta = (TA[0,:,tj,ti])
        loc_dic = (DIC[0,:,tj,ti])
        loc_temp = TEMP[0,:,tj,ti]
        loc_sal = SAL[0,:,tj,ti]

        TA_stor[i,:] = loc_ta
        DIC_stor[i,:] = loc_dic
        TEMP_stor[i,:] = loc_temp
        SAL_stor[i,:] = loc_sal

    ncname = tplcf + '_2015TS.nc'
    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('model_output')
    #g.createDimension('days', len(NO3_mod))
    g.createDimension('days', 365)
    g.createDimension('depths', 40)
    ts = g.createVariable('local_TA','f4',('days','depths'))
    ts[:] = TA_stor
    ts2 = g.createVariable('local_DIC','f4',('days','depths'))
    ts2[:] = DIC_stor
    ts3 = g.createVariable('local_TEMP','f4',('days','depths'))
    ts3[:] = TEMP_stor
    ts4 = g.createVariable('local_SAL','f4',('days','depths'))
    ts4[:] = SAL_stor
    f.close()