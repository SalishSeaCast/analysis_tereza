
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use('seaborn-whitegrid')
import netCDF4 as nc
import cmocean as cm
import glob
import numpy as np
from salishsea_tools import (
    viz_tools,
)
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
import CCCma
import CCCma_stations as cs
from matplotlib import reload
import arrow

def range_analyzer(start,end,varname):
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)
    
    print('days: '+str(len(arrow_array)))
    dayslen = len(arrow_array)
    print('stations: '+str(len(cs.STATIONS)))
    statlen = len(cs.STATIONS)
    
    #DIC,TA,NIT,Sil
    
    storar_mn = np.zeros([statlen,dayslen])
    storar_std = np.zeros([statlen,dayslen])
    
    for i in range(0,len(arrow_array)):
        run_date = arrow_array[i][0]
        ddmmmyy = run_date.format('DDMMMYY').lower()
        humandate = run_date.format('MMM DD, YYYY')
        yyyymmdd = run_date.format('YYYYMMDD')

        print('ANALYZING',humandate)
        #change this if you need to change strings

        carp1 = f'/results2/SalishSea/hindcast.201812_annex/{ddmmmyy}/SalishSea_1d_{yyyymmdd}_{yyyymmdd}_carp_T.nc'
        #grid1 = f'/results2/SalishSea/hindcast.201812_annex/{ddmmmyy}/SalishSea_1d_{yyyymmdd}_{yyyymmdd}_grid_T.nc'
        carp = nc.Dataset(carp1)
        #grid = nc.Dataset(grid1)
        #print(carp)
        var = carp.variables[varname]

        for s in cs.STATIONS:
            #print(s)
            tx = (cs.STATIONS[s]['x'])
            ty = cs.STATIONS[s]['y']
            serno = cs.STATIONS[s]['serialno']
            tdat = var[0,0:5,ty:ty+20,tx:tx+20]
            tdat[tdat==0]=np.nan
            storar_mn[serno,i] = np.nanmean(tdat)
            storar_std[serno,i] = np.nanstd(tdat)
            #print(storar[serno,i])
            
    return storar_mn, storar_std

def range_analyzerptrc(start,end,varname):
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)
    
    print('days: '+str(len(arrow_array)))
    dayslen = len(arrow_array)
    print('stations: '+str(len(cs.STATIONS)))
    statlen = len(cs.STATIONS)
    
    storar_mn = np.zeros([statlen,dayslen])
    storar_std = np.zeros([statlen,dayslen])
    
    for i in range(0,len(arrow_array)):
        run_date = arrow_array[i][0]
        ddmmmyy = run_date.format('DDMMMYY').lower()
        humandate = run_date.format('MMM DD, YYYY')
        yyyymmdd = run_date.format('YYYYMMDD')

        print('ANALYZING',humandate)
        carp1 = f'/results2/SalishSea/hindcast.201812_annex/{ddmmmyy}/SalishSea_1d_{yyyymmdd}_{yyyymmdd}_ptrc_T.nc'
        carp = nc.Dataset(carp1)
        var = carp.variables[varname]

        for s in cs.STATIONS:
            #print(s)
            tx = (cs.STATIONS[s]['x'])
            ty = cs.STATIONS[s]['y']
            serno = cs.STATIONS[s]['serialno']
            tdat = var[0,0:5,ty:ty+20,tx:tx+20]
            tdat[tdat==0]=np.nan
            storar_mn[serno,i] = np.nanmean(tdat)
            storar_std[serno,i] = np.nanstd(tdat)
            #print(storar[serno,i])          
    return storar_mn, storar_std

days = 365
start = '2018-01-01'
end = '2018-12-31'
                                      
print('NO3')
NO3ar, NO3ar_std = range_analyzerptrc(start,end,'nitrate')


ncname = 'NO3_surf5m_050117_080117.nc'
f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('surf_timeseries')
g.createDimension('stn',len(cs.STATIONS))
g.createDimension('time', days)
ts = g.createVariable('timeseries','f4',('stn','time'))
ts[:] = SIar[:]
ts2 = g.createVariable('timeseries_std','f4',('stn','time'))
ts2[:] = SIar_std[:]
f.close()
                                      
print('NH4')
DIATar, DIATar_std = range_analyzerptrc(start,end,'ammonium')


ncname = 'NH4_surf5m_050117_080117.nc'
f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
g = f.createGroup('surf_timeseries')
g.createDimension('stn',len(cs.STATIONS))
g.createDimension('time', days)
ts = g.createVariable('timeseries','f4',('stn','time'))
ts[:] = DIATar[:]
ts2 = g.createVariable('timeseries_std','f4',('stn','time'))
ts2[:] = DIATar_std[:]
f.close()                         
