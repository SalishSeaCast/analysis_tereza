import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    
    us = nc.Dataset(resdir + 'u_1d.nc')
    vs = nc.Dataset(resdir + 'v_1d.nc')
    ugrid = us.variables['vozocrtx']
    vgrid = vs.variables['vomecrty']
    
    thistime = nc.Dataset(resdir + 'timecount_1d.nc')
    timesteps = thistime.variables['time_counter']
    
    return ugrid, vgrid, timesteps, zlevels, depthus

            
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



def velarrow_drawer_DAY(ugrid, vgrid, t, zlevel, dirstr, figtit,depthus,indexer):
# t, zlevel = 0, 0
    bathy = nc.Dataset('/data/dlatorne/MEOPAR/NEMO-forcing/grid/bathy_meter_SalishSea2.nc')
    y_slice = np.arange(180, 500)
    x_slice = np.arange(100, 398)


    arrow_step = 3
    y_slice_a = y_slice[::arrow_step]
    x_slice_a = x_slice[::arrow_step]

    ugrid_tzyx = np.ma.masked_values(ugrid[t, zlevel, y_slice_a, x_slice_a], 0)
    vgrid_tzyx = np.ma.masked_values(vgrid[t, zlevel, y_slice_a, x_slice_a], 0)
    u_tzyx, v_tzyx = viz_tools.unstagger(ugrid_tzyx, vgrid_tzyx)
    speeds = np.sqrt(np.square(u_tzyx) + np.square(v_tzyx))

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    viz_tools.set_aspect(ax)
    quiver = ax.quiver(
        x_slice_a[1:], y_slice_a[1:], u_tzyx, v_tzyx, speeds,
        pivot='mid', cmap=cm.cm.amp, width=0.005)
    viz_tools.plot_land_mask(ax, bathy, xslice=x_slice, yslice=y_slice)

    ax.set_xlim(x_slice[0], x_slice[-1])
    ax.set_ylim(y_slice[0], y_slice[-1])
    ax.grid()

    ax.set_xlabel('x Index')
    ax.set_ylabel('y Index')
    if zlevel == 0:
        tit1 = 'Surface Currents, '
    else:
        tdu = depthus.variables['depthu']
        td = tdu[zlevel]
        tit1 = 'Currents at depth ' + str(td) + ', '
    t_index = indexer + t
    
    
    tit2 = datestring_spitter_DAY(t_index)
    
    si = str(t_index)
    if len(si) == 1:
        lsi = '000' + si
    if len(si) == 2:
        lsi = '00' + si
    if len(si) == 3:
        lsi = '0' + si
    if len(si) == 4:
        lsi = si
    
    tit = tit1 + tit2
    
    ax.set_title(tit, fontsize = 16)
    ax.quiverkey(quiver, 252, 302, 0.5, '0.5 m/s', coordinates='data', color='white', labelcolor='white')
    
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    plt.close(fig)
    
def MRAP(resdir,dirstr,figtit,segment,zlevel):

    day_indexer = [0,15,14,16,15,15,15,15]    
    indexer = (sum(day_indexer[0:segment]))
    print('index')
    print(indexer)
    
    ugrid, vgrid, timesteps, zlevels, depthus  = dat_retrieve_DAY(resdir)
    print('number of frames')
    print(timesteps.shape)
    
    w = (timesteps.shape)
    timer = w[0]
    print(timer)
    for t in range(0,timer):
        velarrow_drawer_DAY(ugrid, vgrid, t, zlevel, dirstr, figtit,depthus,indexer)
        
