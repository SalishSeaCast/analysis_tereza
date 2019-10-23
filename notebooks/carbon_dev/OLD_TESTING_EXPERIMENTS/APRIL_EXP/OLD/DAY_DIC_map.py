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


def nice_42_plot(plotdat, t, v_min, v_max, tcmap, dirstr, figtit, indexer):

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2, figsize = (13, 22))
    land_colour = 'slategrey'
    physdat = nc.Dataset('/data/tjarniko/results/may10_a1/tn_1h.nc')
    zlevels = physdat.variables['deptht']
    
    # Define the component slice to plot
    zmax, ylocn = 41, 423
    section_slice = np.arange(208, 320)
    pdat = np.ma.masked_values(plotdat[t,:,424,section_slice],0)
    cmap = tcmap
    mesh = ax1.pcolormesh(
        section_slice[:], zlevels[:zmax], pdat,
        cmap=cmap, vmin=v_min, vmax=v_max,
    )
    ax1.invert_yaxis()
    cbar = fig.colorbar(mesh, ax=ax1)
    cbar.set_label('DIC μmol/l', fontsize = 15)

    # Axes labels and title
    ax1.set_xlabel('x Index', fontsize = 15)
    ax1.set_ylabel('depth (m)', fontsize = 15)
    #axl.set_title(tit1)

    # Axes limits and grid
    ax1.set_xlim(section_slice[1], section_slice[-1])
    ax1.set_ylim(zlevels[zmax - 2] + 10, 0)
    ax1.set_facecolor(land_colour)
    ax1.grid()
    ax1.set_title('Crossection of Stn. 42', fontsize = 15)

    # Define surface current magnitude slice
    x_slice = np.arange(0, 398)
    y_slice = np.arange(0, 898)
    line_s = np.arange(208, 320)

    surf_dat =  np.ma.masked_values(plotdat[t, 0, y_slice, x_slice], 0)
    
    viz_tools.set_aspect(ax2)
    ax2.plot(
        line_s, 424*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='black',
        label='Section Line',
    )
    
    cmap.set_bad(land_colour)

    mesh = ax2.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)
    ax2.set_xlabel('')
    ax2.set_ylabel('')
    ax2.set_xticks([])
    ax2.set_yticks([])
    #axr.set_title(tit2)
    legend = ax2.legend(loc='best', fancybox=True, framealpha=0.25)
    ax2.grid()
    
    ####
    
    # Define the component slice to plot
    zmax, ylocn = 41, 647
    section_slice = np.arange(120, 220)
    pdat = np.ma.masked_values(plotdat[t,:,ylocn,section_slice],0)
    cmap = tcmap
    mesh = ax3.pcolormesh(
        section_slice[:], zlevels[:zmax], pdat,
        cmap=cmap, vmin=v_min, vmax=v_max,
    )
    ax3.invert_yaxis()
    cbar = fig.colorbar(mesh, ax=ax3)
    cbar.set_label('DIC μmol/l', fontsize = 15)

    # Axes labels and title
    ax3.set_xlabel('x Index', fontsize = 15)
    ax3.set_ylabel('depth (m)', fontsize = 15)
    #axl.set_title(tit1)

    # Axes limits and grid
    ax3.set_xlim(section_slice[1], section_slice[-1])
    ax3.set_ylim(zlevels[zmax - 2] + 10, 0)
    ax3.set_facecolor(land_colour)
    ax3.grid()
    ax3.set_title('Crossection of Stn. 12', fontsize = 15)

    # Define surface current magnitude slice
    x_slice = np.arange(0, 398)
    y_slice = np.arange(0, 898)
    line_s = np.arange(120, 220)

    surf_dat =  np.ma.masked_values(plotdat[t, 0, y_slice, x_slice], 0)
    
    viz_tools.set_aspect(ax4)
    ax4.plot(
        line_s, 647*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='black',
        label='Section Line',
    )
    
    cmap.set_bad(land_colour)

    mesh = ax4.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)
    ax4.set_xlabel('')
    ax4.set_ylabel('')
    ax4.set_xticks([])
    ax4.set_yticks([])
    #axr.set_title(tit2)
    legend = ax4.legend(loc='best', fancybox=True, framealpha=0.25)
    ax4.grid()
    
    
    #####
    
        # Define the component slice to plot
    zmax, ylocn = 41, 505
    section_slice = np.arange(200, 300)
    pdat = np.ma.masked_values(plotdat[t,:,ylocn,section_slice],0)
    cmap = tcmap
    mesh = ax5.pcolormesh(
        section_slice[:], zlevels[:zmax], pdat,
        cmap=cmap, vmin=v_min, vmax=v_max,
    )
    ax5.invert_yaxis()
    cbar = fig.colorbar(mesh, ax=ax5)
    cbar.set_label('DIC μmol/l', fontsize = 15)

    # Axes labels and title
    ax5.set_xlabel('x Index', fontsize = 15)
    ax5.set_ylabel('depth (m)', fontsize = 15)
    #axl.set_title(tit1)

    # Axes limits and grid
    ax5.set_xlim(section_slice[1], section_slice[-1])
    ax5.set_ylim(zlevels[zmax -2] + 20, 0)
    ax5.set_facecolor(land_colour)
    ax5.grid()
    ax5.set_title('Crossection of Stn. 27', fontsize = 15)

    # Define surface current magnitude slice
    x_slice = np.arange(0, 398)
    y_slice = np.arange(0, 898)
    line_s = np.arange(200, 300)

    surf_dat =  np.ma.masked_values(plotdat[t, 0, y_slice, x_slice], 0)
    
    viz_tools.set_aspect(ax6)
    ax6.plot(
        line_s, ylocn*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='black',
        label='Section Line',
    )
    
    cmap.set_bad(land_colour)

    mesh = ax6.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)
    ax6.set_xlabel('')
    ax6.set_ylabel('')
    ax6.set_xticks([])
    ax6.set_yticks([])
    #axr.set_title(tit2)
    legend = ax6.legend(loc='best', fancybox=True, framealpha=0.25)
    ax6.grid()
    
 
    

    t_index = indexer + t
    si = str(t_index)
    if len(si) == 1:
        lsi = '000' + si
    if len(si) == 2:
        lsi = '00' + si
    if len(si) == 3:
        lsi = '0' + si
    if len(si) == 4:
        lsi = si
    
    
    
    
    tit2 = datestring_spitter_DAY(t_index)
    
    tit = 'DIC ' + tit2
    plt.suptitle(tit, fontsize = 20)
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    #plt.show()
    plt.close(fig)

def nice_42_plot_orig(plotdat, t, v_min, v_max, tcmap, dirstr, figtit, indexer):

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(15, 8 ))
    land_colour = 'slategrey'

    physdat = nc.Dataset('/data/tjarniko/results/may10_a1/tn_1h.nc')
    zlevels = physdat.variables['deptht']

    # Define the component slice to plot
    zmax, ylocn = 41, 424
    section_slice = np.arange(208, 293)

    pdat = np.ma.masked_values(plotdat[t,:,424,section_slice],0)
    
    cmap = tcmap
    #cmap.set_bad(land_colour)
    #cmap.set_bad('whitesmoke')
    mesh = axl.pcolormesh(
        section_slice[:], zlevels[:zmax], pdat,
        cmap=cmap, vmin=v_min, vmax=v_max,
    )
    axl.invert_yaxis()
    cbar = fig.colorbar(mesh, ax=axl)
    cbar.set_label('DIC μmol/l', fontsize = 15)

    # Axes labels and title
    axl.set_xlabel('x Index', fontsize = 15)
    axl.set_ylabel('depth (m)', fontsize = 15)
    #axl.set_title(tit1)

    # Axes limits and grid
    axl.set_xlim(section_slice[1], section_slice[-1])
    axl.set_ylim(zlevels[zmax - 2] + 10, 0)
    axl.set_facecolor(land_colour)
    axl.grid()

    # Define surface current magnitude slice
    x_slice = np.arange(0, 398)
    y_slice = np.arange(0, 898)
    line_s = np.arange(0,398)

    
    surf_dat =  np.ma.masked_values(plotdat[t, 0, y_slice, x_slice], 0)
    
    viz_tools.set_aspect(axr)
    axr.plot(
        line_s, 424*np.ones_like(line_s),
        linestyle='solid', linewidth=3, color='black',
        label='Section Line',
    )
    
    cmap.set_bad(land_colour)

    mesh = axr.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)

    axr.set_xlabel('')
    axr.set_ylabel('')
    axr.set_xticks([])
    axr.set_yticks([])
    #axr.set_title(tit2)
    legend = axr.legend(loc='best', fancybox=True, framealpha=0.25)
    axr.grid()

    t_index = indexer + t
    si = str(t_index)
    if len(si) == 1:
        lsi = '000' + si
    if len(si) == 2:
        lsi = '00' + si
    if len(si) == 3:
        lsi = '0' + si
    if len(si) == 4:
        lsi = si
    
    
    
    
    tit2 = datestring_spitter_DAY(t_index)
    
    tit = 'DIC ' + tit2
    plt.suptitle(tit, fontsize = 20)
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    plt.close(fig)
    
def MRAP(resdir,dirstr,figtit,segment):

    day_indexer = [0,15,15,14,16,15,15,15]   
    indexer = (sum(day_indexer[0:segment]))
    print('index')
    print(indexer)
    
    trac, ptrac, timesteps  = dat_retrieve_DAY(resdir)
    print('number of frames')
    print(timesteps.shape)
    
    dic = trac.variables['dissolved_inorganic_carbon']
    w = (timesteps.shape)
    
    tcmap = cm.cm.matter
    v_min = 1800
    v_max = 2200
    timer = w[0]
    print(timer)
    for t in range(0,timer):
        nice_42_plot(dic, t, v_min, v_max, tcmap, dirstr, figtit, indexer)
        

