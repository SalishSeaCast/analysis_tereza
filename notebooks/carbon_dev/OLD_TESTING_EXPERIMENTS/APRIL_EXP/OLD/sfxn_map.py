import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    
    us = nc.Dataset(resdir + 'u_1h.nc')
    vs = nc.Dataset(resdir + 'v_1h.nc')
    ugrid = us.variables['vozocrtx']
    vgrid = vs.variables['vomecrty']
    
    thistime = nc.Dataset(resdir + 'timecount_1h.nc')
    timesteps = thistime.variables['time_counter']
    
    return ugrid, vgrid, timesteps, zlevels, depthus
        
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


def streamfxn_drawer(ugrid, vgrid, t, zlevel, dirstr, figtit, depthus, indexer):
    bathy = nc.Dataset('/data/dlatorne/MEOPAR/NEMO-forcing/grid/bathy_meter_SalishSea2.nc')
    # t, zlevel = 0, 0
    y_slice = np.arange(180, 500)
    x_slice = np.arange(100, 398)

    ugrid_tzyx = np.ma.masked_values(ugrid[t, zlevel, y_slice, x_slice], 0)
    vgrid_tzyx = np.ma.masked_values(vgrid[t, zlevel, y_slice, x_slice], 0)
    u_tzyx, v_tzyx = viz_tools.unstagger(ugrid_tzyx, vgrid_tzyx)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    density = 10
    # for ax, density in zip(axs, densities):
    viz_tools.set_aspect(ax)

    ax.streamplot(
        x_slice[1:], y_slice[1:], u_tzyx, v_tzyx,
        density=density,
    )
    viz_tools.plot_land_mask(ax, bathy, xslice=x_slice, yslice=y_slice)

    ax.set_xlim(x_slice[0], x_slice[-1])
    ax.set_ylim(y_slice[0], y_slice[-1])
    ax.grid()

    ax.set_xlabel('x Index')
    ax.set_ylabel('y Index')
    if zlevel == 0:
        tit1 = 'Surface Streamlines, '
    else:
        tdu = depthus.variables['depthu']
        td = tdu[zlevel]
        tit1 = 'Streamlines at depth ' + str(td) + ', '
    t_index = indexer + t
    
    tit2 = datestring_spitter(t_index)
    
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

    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    plt.close(fig)


    
def MRAP(resdir,dirstr,figtit,segment,zlevel):

    segment_lister = [0,360,720,1056,1440,1800,2160,2520]
    indexer = (segment_lister[segment-1])
    print('index')
    print(indexer)
    ugrid, vgrid, timesteps, zlevels, depthus = dat_retrieve(resdir)
    print('number of frames')
    print(timesteps.shape)
    
    w = (timesteps.shape)
    timer = w[0]
    
    for t in range(0,timer):
        streamfxn_drawer(ugrid, vgrid, t, zlevel, dirstr, figtit,depthus,indexer)

    
 

