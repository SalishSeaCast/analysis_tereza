import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn')
from salishsea_tools import visualisations as vis
from salishsea_tools import (teos_tools, tidetools, viz_tools)
import cmocean as cm
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

def dat_retrieve(resdir):
    
    resdir = '/data/tjarniko/results/' + resdir
    
    print(resdir)
    
    depthus = nc.Dataset('/data/tjarniko/results/may10_a1/zlevels_1h.nc')
    zlevels = depthus.variables['depthu']
    
    DIC =  nc.Dataset('/data/tjarniko/results/may10_a1/DIC_1h.nc')
    OXY =  nc.Dataset('/data/tjarniko/results/may10_a1/OXY_1h.nc')
    sn =  nc.Dataset('/data/tjarniko/results/may10_a1/sn_1h.nc')


    thistime = nc.Dataset(resdir + 'timecount_1h.nc')
    timesteps = thistime.variables['time_counter']
    
    return DIC, OXY, sn, thistime, timesteps
        
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


def dicoxysal_prof(t,DIC,OXY,sn,dirstr,figtit,indexer):

    #specific to station 42
    t_y = 426
    t_x = 255
    dend = 33
    
    
    sal_prof = sn.variables['vosaline'][0,:,t_y,t_x]
    DIC_prof = DIC.variables['dissolved_inorganic_carbon'][0,:,t_y,t_x]
    OXY_prof = OXY.variables['dissolved_oxygen'][0,:,t_y,t_x]
    depth = sn.variables['deptht'][:]
    # Three-panel plot
    fig, (ax2, ax3, ax4) = plt.subplots(figsize=(14.0, 9.0) , nrows=1, ncols=3, sharey=True)
    # Temperature
    ax2.plot(sal_prof[0:dend],depth[0:dend],'o-')
    ax2.set_ylabel('Depth (m)', fontsize = 14)
    ax2.set_ylim([0,255])
    ax2.set_ylim(ax2.get_ylim()[::-1]) #this reverses the yaxis (i.e. deep at the bottom)
    ax2.set_xlabel('Salinity', fontsize = 14)
    #ax2.xaxis.set_label_position('top') # this moves the label to the top
    #ax2.xaxis.set_ticks_position('top') # this moves the ticks to the top
    # Salinity
    ax3.plot(DIC_prof[0:dend],depth[0:dend],'o-r')
    ax3.set_xlabel('DIC μmol/kg', fontsize = 14)
    #ax3.xaxis.set_label_position('top') # this moves the label to the top
    #ax3.xaxis.set_ticks_position('top') # this moves the ticks to the top
    ax3.yaxis.set_visible(True) # This erases the y ticks
    # Fluorescence
    ax4.plot(OXY_prof[0:dend],depth[0:dend],'o-g')
    ax4.set_xlabel('oxygen ml/l', fontsize = 14)
    #ax4.xaxis.set_label_position('top') # this moves the label to the top
    #ax4.xaxis.set_ticks_position('top') # this moves the ticks to the top
    ax4.yaxis.set_visible(True) # This erases the y ticks
    tit1 = 'Depth Prof., Stn 42. '
    t_index = indexer + t
    tit2 = datestring_spitter(t_index)
    fig.suptitle(tit1 + tit2, fontsize = 18)
    
    
    si = str(t_index)
    if len(si) == 1:
        lsi = '000' + si
    if len(si) == 2:
        lsi = '00' + si
    if len(si) == 3:
        lsi = '0' + si
    if len(si) == 4:
        lsi = si
        
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    plt.close(fig)

    
def MRAP(resdir,dirstr,figtit,segment):

    segment_lister = [0,360,720,1056,1440,1800,2160,2520]
    indexer = (segment_lister[segment-1])
    print('index')
    print(indexer)
    DIC, OXY, sn, thistime, timesteps = dat_retrieve(resdir)
    print('number of frames')
    print(timesteps.shape)
    
    w = (timesteps.shape)
    timer = w[0]
    
    for t in range(0,timer):
        dicoxysal_prof(t,DIC,OXY,sn,dirstr,figtit,indexer)
    
 

