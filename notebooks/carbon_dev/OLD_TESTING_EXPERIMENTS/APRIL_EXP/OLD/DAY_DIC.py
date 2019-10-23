#sample use
#resdir = 'may10_a1/'
#zlevel = 0
#segment = 1
#dirstr = './APR_dp_day/'
#figtit = 'DP_day_'
#Ddp.MRAP_day(resdir,dirstr,figtit,segment)

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
import glob 

def dat_retrieve_DAY(resdir):
    
    resdir = '/data/tjarniko/results/' + resdir
    
    print(resdir)
    
    depthus = nc.Dataset('/data/tjarniko/results/may10_a1/zlevels_1h.nc')
    zlevels = depthus.variables['depthu']
    w = glob.glob(resdir + 'SKOG_1d*ptrc_T.nc')
    w2 = glob.glob(resdir + 'SKOG_1d*grid_T.nc')
    trac =  nc.Dataset(w[0])
    ptrac =  nc.Dataset(w2[0])                  
    timesteps = trac.variables['time_counter']
    
    return trac, ptrac, timesteps
        
def datestring_spitter_DAY(days_since):
    days = [31,29,31,30]
    
    month = 'walrus'
    if (days_since < (sum(days[0:4]))):
        month = 'April'
        left = (days_since - sum(days[0:3]))
        day = np.floor(left)
        day = day + 1
        hour = left%24
    if (days_since < (sum(days[0:3]))):
        month = 'March'
        left = (days_since - sum(days[0:2]))
        day = np.floor(left)
        day = day + 1
    
    if (days_since < (sum(days[0:2]))):
        month = 'February'
        left = (days_since - sum(days[0:1]))
        day = np.floor(left)
        
        day = day + 1
        hour = left%24
    if (days_since < (sum(days[0:1]))):
        month = 'January'
        left = (days_since - sum(days[0:0]))
        day = np.floor(left)
        day = day + 1

    datestr = month + ' ' + str(int(day)) + ', 2016, (' + str(days_since) + ' days since January 1, 2016)'
    print(datestr)
    return datestr


def dicoxysal_prof(t,trac,ptrac,dirstr,figtit,indexer):

    #specific to station 42
    t_y = 426
    t_x = 255
    dend = 33
    
    print('t is')
    print(t)
    sal_prof = ptrac.variables['vosaline'][t,:,t_y,t_x]
    DIC_prof = trac.variables['dissolved_inorganic_carbon'][t,:,t_y,t_x]
    OXY_prof = trac.variables['dissolved_oxygen'][t,:,t_y,t_x]
    depth = ptrac.variables['deptht'][:]
    # Three-panel plot
    fig, (ax2, ax3, ax4) = plt.subplots(figsize=(14.0, 9.0) , nrows=1, ncols=3, sharey=True)
    # Temperature
    ax2.plot(sal_prof[0:dend],depth[0:dend],'o-')
    ax2.set_ylabel('Depth (m)', fontsize = 14)
    ax2.set_ylim([0,255])
    ax2.set_ylim(ax2.get_ylim()[::-1]) #this reverses the yaxis (i.e. deep at the bottom)
    ax2.set_xlim([26,32])
    ax2.set_xlabel('Salinity', fontsize = 14)
    #ax2.xaxis.set_label_position('top') # this moves the label to the top
    #ax2.xaxis.set_ticks_position('top') # this moves the ticks to the top
    # Salinity
    ax3.plot(DIC_prof[0:dend],depth[0:dend],'o-r')
    ax3.set_xlim([1800,2200])
    ax3.set_xlabel('DIC μmol/kg', fontsize = 14)
    #ax3.xaxis.set_label_position('top') # this moves the label to the top
    #ax3.xaxis.set_ticks_position('top') # this moves the ticks to the top
    ax3.yaxis.set_visible(True) # This erases the y ticks
    # Fluorescence
    ax4.plot(OXY_prof[0:dend],depth[0:dend],'o-g')
    ax4.set_xlim([3,10])
    ax4.set_xlabel('oxygen ml/l', fontsize = 14)
    #ax4.xaxis.set_label_position('top') # this moves the label to the top
    #ax4.xaxis.set_ticks_position('top') # this moves the ticks to the top
    ax4.yaxis.set_visible(True) # This erases the y ticks
    tit1 = 'Depth Prof., Stn 42. '
    t_index = indexer + t
    tit2 = datestring_spitter_DAY(t_index)
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

    
def MRAP_day(resdir,dirstr,figtit,segment):

    day_indexer = [0,15,15,14,16,15,15,15]    
    indexer = (sum(day_indexer[0:segment])) 
    print('index')
    print(indexer)
    
    trac, ptrac, timesteps = dat_retrieve_DAY(resdir)
    print('number of frames')
    print(timesteps.shape)
    
    w = (timesteps.shape)
    timer = w[0]
    
    for w in range(0,timer):
        print('no look I pass it the right effing t')
        print(w)
        dicoxysal_prof(w,trac,ptrac,dirstr,figtit,indexer)
    
 


