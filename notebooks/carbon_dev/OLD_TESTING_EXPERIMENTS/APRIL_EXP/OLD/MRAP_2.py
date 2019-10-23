import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from salishsea_tools import visualisations as vis
import cmocean as cm
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

def dat_retrieve(resdir):
    
    resdir = '/data/tjarniko/results/' + resdir
    
    print(resdir)

    DIC = nc.Dataset(resdir + 'DIC_1h.nc')
    TA = nc.Dataset(resdir + 'TA_1h.nc')
    sn = nc.Dataset(resdir + 'sn_1h.nc')
    tn = nc.Dataset(resdir + 'tn_1h.nc')
    OXY = nc.Dataset(resdir + 'OXY_1h.nc')
    
    DIC_dat = DIC.variables['dissolved_inorganic_carbon']
    TA_dat = TA.variables['total_alkalinity']
    sn_dat = sn.variables['vosaline']
    tn_dat = tn.variables['votemper']
    OXY_dat = OXY.variables['dissolved_oxygen']

    print(sn_dat.shape)
    no_frames = sn_dat.shape[0]
    return no_frames, DIC_dat, TA_dat, sn_dat, tn_dat, OXY_dat
    
def thalweg_animator(plotdat,tstart,tend,vmin,vmax,stepsize,tit,figtit,dirstr,indexer,t_cmap,clabel):
    "TESTED"

    
    bathy = nc.Dataset('/data/tjarniko/MEOPAR/grid/bathymetry_201702.nc')
    mesh = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
    
    for i in range(tstart,tend):
        td = plotdat[i,:,:,:]
        cmap = t_cmap
        fig,ax = plt.subplots(1,1,figsize=(15,5))

        cbar = vis.contour_thalweg(ax, td, bathy, mesh, np.arange(vmin, vmax, stepsize), cmap = t_cmap)
        cbar.set_label(clabel, fontsize = 20)
        t_i = indexer + i
        si = str(t_i)
        #print(len(si))
        if len(si) == 1:
            lsi = '000' + si
        if len(si) == 2:
            lsi = '00' + si
        if len(si) == 3:
            lsi = '0' + si
        if len(si) == 4:
            lsi = si
        print(lsi)
        fig.suptitle(datestring_spitter(t_i))
        total_fig = dirstr+figtit+str(lsi) + '.png'
        fig.savefig(total_fig)
        plt.close(fig)
        
def datestring_spitter(hrs_since):
    days = [31,29,31,30]
    h=24
    dayhr = h*np.array(days)
    #print(dayhr)

    month = 'walrus'
    day = 'walrus'
    hour = 'walrus'
    
    if (hrs_since < (sum(dayhr[0:4]))):
        month = 'April'
        left = (hrs_since - sum(dayhr[0:3]))
        day = np.floor(left/24)
        day = day + 1
        hour = left%24
    if (hrs_since < (sum(dayhr[0:3]))):
        month = 'March'
        left = (hrs_since - sum(dayhr[0:2]))
        day = np.floor(left/24)
        day = day + 1
        hour = left%24
    if (hrs_since < (sum(dayhr[0:2]))):
        month = 'February'
        left = (hrs_since - sum(dayhr[0:1]))
        day = np.floor(left/24)
        day = day + 1
        hour = left%24
    if (hrs_since < (sum(dayhr[0:1]))):
        month = 'January'
        left = hrs_since 
        day = np.floor(left/24)
        day = day + 1
        hour = left%24

    hour = str(hour)
    if len(hour) == 1:
        hour = '0' + hour + ':00'
    else: 
        hour = hour + ':00'
    
    day = str(int(day))
    
    datestr = month + ' ' + day + ', 2016, ' + hour + '  (' + str(hrs_since) + ' hours since January 1, 2016)'
    print(datestr)
    return datestr

    
def MRAP(resdir,DIC,TA,OXY,sn,tn,dirstr,segment):

    segment_lister = [0,360,720,1056,1440,1800,2160,2520]
    indexer = (segment_lister[segment-1])
    print('index')
    print(indexer)
    no_frames, DIC_dat, TA_dat, sn_dat, tn_dat, OXY_dat = dat_retrieve(resdir)
    print('number of frames')
    print(no_frames)
#     no_frames = 2
# #     print('hour we are at')
# #     print(indexer)
    if (DIC == 'DICy'):
        print('DIC plot')
        t_cmap = cm.cm.matter
        vmin = 1800
        vmax = 2400
        stepsize = 12
        tit = 'DIC along thalweg, '
        tit2 = ' hours since January 1, 2016'
        figtit = 'DIC_1hr_'
        clabel = 'DIC, μmol/kg'
        thalweg_animator(DIC_dat,0,no_frames,vmin,vmax,stepsize,tit,figtit,dirstr,indexer,t_cmap,clabel)
        
    if (TA == 'TAy'):
        print('TA plot')
        t_cmap = cm.cm.matter
        vmin = 1800
        vmax = 2400
        stepsize = 12
        tit = 'TA along thalweg, '
        tit2 = ' hours since January 1, 2016'
        figtit = 'TA_1hr_'
        clabel = 'TA, μmol/kg'
        thalweg_animator(TA_dat,0,no_frames,vmin,vmax,stepsize,tit,figtit,dirstr,indexer,t_cmap,clabel)       
               
    if (OXY == 'OXYy'):
        print('OXY plot')
        t_cmap = cm.cm.dense
        vmin = 0
        vmax = 12
        stepsize = 0.1
        tit = 'Oxygen along thalweg, '
        tit2 = ' hours since January 1, 2016'
        figtit = 'OXY_1hr_'
        clabel = 'oxygen, μmol/kg'
        thalweg_animator(OXY_dat,0,no_frames,vmin,vmax,stepsize,tit,figtit,dirstr,indexer,t_cmap,clabel)
        

    if (sn == 'sny'):
        print('sn plot')
        t_cmap = cm.cm.haline
        vmin = 20
        vmax = 35
        stepsize = 0.1
        tit = 'Salinity along thalweg, '
        tit2 = ' hours since January 1, 2016'
        figtit = 'sn_1hr_'
        clabel = 'salinity, psu'
        thalweg_animator(sn_dat,0,no_frames,vmin,vmax,stepsize,tit,figtit,dirstr,indexer,t_cmap,clabel)
        
    if (tn == 'tny'):
        print('tn plot')
        t_cmap = cm.cm.thermal
        vmin = 0
        vmax = 12
        stepsize = 0.1
        tit = 'Temperature along thalweg, '
        tit2 = ' hours since January 1, 2016'
        figtit = 'tn_1hr_'
        clabel = 'temperature (°C)'
        thalweg_animator(tn_dat,0,no_frames,vmin,vmax,stepsize,tit,figtit,dirstr,indexer,t_cmap,clabel)

MRAP('may10_a2/','DICy','TAy','OXYy','sny','tny','./APR_thwg/',2)
