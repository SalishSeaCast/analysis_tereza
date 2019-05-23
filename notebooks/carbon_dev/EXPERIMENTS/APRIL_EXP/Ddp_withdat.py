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

def dat_retrieve_fielddat():
    infil = np.loadtxt('/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/DATASETS/2016_05.txt')
    crid = infil[:,0]
    ln = infil[:,2]
    stn = infil[:,3]
    lat = infil[:,6]
    lon = infil[:,7]
    P = infil[:,8]
    T = infil[:,9]
    S = infil[:,10]
    ox = infil[:,11]
    dic = infil[:,13]
    st_42 = (stn == 42)
    st_12 = (stn == 12)
    st_27 = (stn == 27)
    S = teos_tools.psu_teos(S)
    
    return stn, lat, lon, P, T, S, ox, dic, st_42, st_12, st_27

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


def dicoxysal_prof(t,trac,ptrac,dirstr,figtit,indexer,stn, lat, lon, P, T, S, ox, dic, st_42, st_12, st_27):

    #specific to station 42 - south, red
    y42 = 423
    x42 = 263
    
    #specific to stations 12 and 27 - north, blue
    y12 = 647
    x12 = 168
    
    y27 = 505
    x27 = 246
    
    dend = 33
    
    print('t is')
    print(t)
    sal_prof_42 = ptrac.variables['vosaline'][t,:,y42,x42]
    DIC_prof_42 = trac.variables['dissolved_inorganic_carbon'][t,:,y42,x42]
    OXY_prof_42 = trac.variables['dissolved_oxygen'][t,:,y42,x42]
    
    sal_prof_12 = ptrac.variables['vosaline'][t,:,y12,x12]
    DIC_prof_12 = trac.variables['dissolved_inorganic_carbon'][t,:,y12,x12]
    OXY_prof_12 = trac.variables['dissolved_oxygen'][t,:,y12,x12]
    
    sal_prof_27 = ptrac.variables['vosaline'][t,:,y27,x27]
    DIC_prof_27 = trac.variables['dissolved_inorganic_carbon'][t,:,y27,x27]
    OXY_prof_27 = trac.variables['dissolved_oxygen'][t,:,y27,x27]
    
    depth = ptrac.variables['deptht'][:]
    # Three-panel plot
    fig, (ax2, ax3) = plt.subplots(figsize=(12.0, 9.0) , nrows=1, ncols=2, sharey=True)
    # Temperature
    ax2.plot(S[st_42],P[st_42],linestyle='', marker ='o', color = 'xkcd:red', markersize = 9)    
    ax2.plot(S[st_12],P[st_12],linestyle='', marker ='o', color = 'xkcd:bright blue', markersize = 9)
    ax2.plot(S[st_27],P[st_27],linestyle='', marker ='o', color = 'xkcd:royal blue', markersize = 9)
    ax2.plot(sal_prof_42[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='red')
    ax2.plot(sal_prof_12[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:bright blue')
    ax2.plot(sal_prof_27[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:royal blue')



    ax2.set_ylabel('Depth (m)', fontsize = 18)
    ax2.set_ylim([0,255])
    ax2.set_ylim(ax2.get_ylim()[::-1]) #this reverses the yaxis (i.e. deep at the bottom)
    ax2.set_xlim([26,32])
    ax2.set_xlabel('Salinity', fontsize = 18)
    #ax2.xaxis.set_label_position('top') # this moves the label to the top
    #ax2.xaxis.set_ticks_position('top') # this moves the ticks to the top

    ax3.plot(dic[st_42],P[st_42], linestyle='', marker ='o', color = 'xkcd:red', markersize = 9)
    ax3.plot(dic[st_12],P[st_12], linestyle='', marker ='o', color = 'xkcd:bright blue', markersize = 9)
    ax3.plot(dic[st_27],P[st_27], linestyle='', marker ='o', color = 'xkcd:royal blue', markersize = 9)    
    ax3.plot(DIC_prof_42[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='red')
    ax3.plot(DIC_prof_12[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:bright blue')
    ax3.plot(DIC_prof_27[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:royal blue')



    
    ax3.set_xlim([1800,2200])
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax3.tick_params(axis='both', which='major', labelsize=15)
    ax3.set_xlabel('DIC μmol/kg', fontsize = 18)
    #ax3.xaxis.set_label_position('top') # this moves the label to the top
    #ax3.xaxis.set_ticks_position('top') # this moves the ticks to the top
    ax3.yaxis.set_visible(True) # This erases the y ticks

    ax3.legend(['stn 42','stn 12','stn 27'], fontsize = 18)

    tit1 = 'Depth Profiles of Stns 12, 27, 42, '
    t_index = indexer + t
    tit2 = datestring_spitter_DAY(t_index)
    fig.suptitle(tit1 + tit2, fontsize = 22)
    
    
    si = str(t_index)
    if len(si) == 1:
        lsi = '000' + si
    if len(si) == 2:
        lsi = '00' + si
    if len(si) == 3:
        lsi = '0' + si
    if len(si) == 4:
        lsi = si

    ax3.legend(['SOUTH','NORTH','NORTH-CENTRAL'], fontsize = 18)
    ax3.xaxis.set_tick_params(labelsize=18)
    ax3.yaxis.set_tick_params(labelsize=18)
    ax2.xaxis.set_tick_params(labelsize=18)
    ax2.yaxis.set_tick_params(labelsize=18)    
        
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    plt.close(fig)

    
def dicoxysal_profss(t,trac,ptrac,dirstr,figtit,indexer,stn, lat, lon, P, T, S, ox, dic, st_42, st_12, st_27):

    #specific to station 42 - south, red
    y42 = 423
    x42 = 263
    
    #specific to stations 12 and 27 - north, blue
    y12 = 647
    x12 = 168
    
    y27 = 505
    x27 = 246
    
    dend = 33
    
    print('t is')
    print(t)
    sal_prof_42 = ptrac.variables['vosaline'][t,:,y42,x42]
    DIC_prof_42 = trac.variables['dissolved_inorganic_carbon'][t,:,y42,x42]
    OXY_prof_42 = trac.variables['dissolved_oxygen'][t,:,y42,x42]
    
    sal_prof_12 = ptrac.variables['vosaline'][t,:,y12,x12]
    DIC_prof_12 = trac.variables['dissolved_inorganic_carbon'][t,:,y12,x12]
    OXY_prof_12 = trac.variables['dissolved_oxygen'][t,:,y12,x12]
    
    sal_prof_27 = ptrac.variables['vosaline'][t,:,y27,x27]
    DIC_prof_27 = trac.variables['dissolved_inorganic_carbon'][t,:,y27,x27]
    OXY_prof_27 = trac.variables['dissolved_oxygen'][t,:,y27,x27]
    
    depth = ptrac.variables['deptht'][:]
    # Three-panel plot
    fig, (ax2, ax3) = plt.subplots(figsize=(12.0, 9.0) , nrows=1, ncols=2, sharey=True)
    # Temperature
    ax2.plot(S[st_42],P[st_42],linestyle='', marker ='o', color = 'xkcd:red', markersize = 11)    
    ax2.plot(S[st_12],P[st_12],linestyle='', marker ='o', color = 'xkcd:bright blue', markersize = 11)
    ax2.plot(S[st_27],P[st_27],linestyle='', marker ='o', color = 'xkcd:royal blue', markersize = 11)
    ax2.plot(sal_prof_42[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='red')
    ax2.plot(sal_prof_12[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:bright blue')
    ax2.plot(sal_prof_27[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:royal blue')



    ax2.set_ylabel('Depth (m)', fontsize = 18)
    ax2.set_ylim([20,150])
    ax2.set_ylim(ax2.get_ylim()[::-1]) #this reverses the yaxis (i.e. deep at the bottom)
    ax2.set_xlim([26,32])
    ax2.set_xlabel('Salinity', fontsize = 18)
    #ax2.xaxis.set_label_position('top') # this moves the label to the top
    #ax2.xaxis.set_ticks_position('top') # this moves the ticks to the top

    ax3.plot(dic[st_42],P[st_42], linestyle='', marker ='o', color = 'xkcd:red', markersize = 11)
    ax3.plot(dic[st_12],P[st_12], linestyle='', marker ='o', color = 'xkcd:bright blue', markersize = 11)
    ax3.plot(dic[st_27],P[st_27], linestyle='', marker ='o', color = 'xkcd:royal blue', markersize = 11)    
    ax3.plot(DIC_prof_42[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='red')
    ax3.plot(DIC_prof_12[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:bright blue')
    ax3.plot(DIC_prof_27[0:dend],depth[0:dend],linestyle='-', linewidth = 3, color='xkcd:royal blue')



    
    ax3.set_xlim([1920,2060])
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax3.tick_params(axis='both', which='major', labelsize=15)
    ax3.set_xlabel('DIC μmol/kg', fontsize = 18)
    #ax3.xaxis.set_label_position('top') # this moves the label to the top
    #ax3.xaxis.set_ticks_position('top') # this moves the ticks to the top
    ax3.yaxis.set_visible(True) # This erases the y ticks

    ax3.legend(['stn 42','stn 12','stn 27'], fontsize = 18)

    tit1 = 'Depth Profiles, '
    t_index = indexer + t
    tit2 = datestring_spitter_DAY(t_index)
    fig.suptitle(tit1 + tit2, fontsize = 22)
    
    
    si = str(t_index)
    if len(si) == 1:
        lsi = '000' + si
    if len(si) == 2:
        lsi = '00' + si
    if len(si) == 3:
        lsi = '0' + si
    if len(si) == 4:
        lsi = si
        
    ax3.legend(['SOUTH','NORTH','NORTH-CENTRAL'], fontsize = 18)
    ax3.xaxis.set_tick_params(labelsize=18)
    ax3.yaxis.set_tick_params(labelsize=18)
    ax2.xaxis.set_tick_params(labelsize=18)
    ax2.yaxis.set_tick_params(labelsize=18)
    
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    plt.close(fig)

    
def MRAP_day(resdir,dirstr,figtit,segment):
    
    stn, lat, lon, P, T, S, ox, dic, st_42, st_12, st_27 = dat_retrieve_fielddat()

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
        dicoxysal_profss(w,trac,ptrac,dirstr,figtit,indexer, stn, lat, lon, P, T, S, ox, dic, st_42, st_12, st_27)
    
 


