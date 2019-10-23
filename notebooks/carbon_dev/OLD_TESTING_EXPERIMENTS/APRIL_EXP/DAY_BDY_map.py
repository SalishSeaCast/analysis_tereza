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
    w3 = glob.glob(resdir + 'SKOG_1d*grid_V.nc')
    trac =  nc.Dataset(w[0])
    ptrac =  nc.Dataset(w2[0])
    vtrac = nc.Dataset(w3[0])
    timesteps = trac.variables['time_counter']
    ratrac = nc.Dataset('/data/tjarniko/results/apr_nces/SKOG_ncraJAN.nc')
    
    return trac, ptrac, vtrac, ratrac, timesteps

            
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

def surf_plot(plotdat, t, v_min, v_max, tcmap, dirstr, figtit, indexer,surtit, deep, cblabel):

    fig, ((ax4)) = plt.subplots(1,1, figsize = (8.5, 12.5))
    land_colour = 'silver'
    physdat = nc.Dataset('/data/tjarniko/results/may10_a1/tn_1h.nc')
    zlevels = physdat.variables['deptht']
    print('depth')
    print(zlevels[deep])
    td = (int(zlevels[deep]))
    print(td)
    cmap = tcmap
    
    ### ax4 is the map
    
    x_slice = np.arange(0, 398)
    y_slice = np.arange(0, 898)
    line_s = np.arange(120, 220)

    surf_dat =  np.ma.masked_values(plotdat[t, deep, :, :], 0)
    #surf_dat =  plotdat[t, deep, y_slice, x_slice]
    
    viz_tools.set_aspect(ax4)
    
    line_s = np.arange(200, 320)
    ax4.plot(
        line_s, 423*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='red',
        label='Section Line',
    )
    i = 263
    j = 423
    ax4.plot(i,j,marker='o',color='xkcd:red',markersize=12, markeredgecolor = 'w', markeredgewidth = 2.0)
    ax4.text(i-150,j,'SOUTH',fontsize = 13)
    
    line_s = np.arange(190, 310)
    ax4.plot(
        line_s, 505*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='xkcd:royal blue',
        label='Section Line',
    )
    i = 246
    j = 505
    ax4.plot(i,j,marker='o',color='xkcd:royal blue',markersize=12, markeredgecolor = 'w', markeredgewidth = 2.0)
    ax4.text(i-250,j,'N.-CENTRAL',fontsize = 13)
    
    
    
    line_s = np.arange(110, 230)
    ax4.plot(
        line_s, 647*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='xkcd:bright blue',
        label='Section Line',
    )
    i = 168
    j = 647
    ax4.plot(i,j,marker='o',color='xkcd:bright blue',markersize=12, markeredgecolor = 'w', markeredgewidth = 2.0)
    print('text?')
    ax4.text(i-150,j,'NORTH',fontsize = 13)
    
    
    cmap.set_bad(land_colour)
    mesh = ax4.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)
    ax4.set_xlabel('')
    ax4.set_ylabel('')
    #ax4.set_xticks([])
    #ax4.set_yticks([])
    #legend = ax4.legend(loc='best', fancybox=True, framealpha=0.25)
    ax4.grid()
    cbar = fig.colorbar(mesh, ax=ax4)
    cbar.set_label(cblabel, fontsize = 15)
    #####
    
    if td < 1:
        ax4.set_title('Surface ' + surtit, fontsize = 20, color = 'black')
    else:
        ax4.set_title(surtit +'at depth ' + str(td) +' m', fontsize = 20, color = 'black')


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
    
    tit = surtit + tit2
    plt.suptitle(tit, fontsize = 20)
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    #plt.show()
    plt.close(fig)
    
def a3stn_plot(plotdat, t, v_min, v_max, tcmap, dirstr, figtit, indexer, surtit, deep, cblabel):

    fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(1,4, figsize = (20, 7))
    land_colour = 'silver'
    physdat = nc.Dataset('/data/tjarniko/results/may10_a1/tn_1h.nc')

    zlevels = physdat.variables['deptht']
    #print('depth')
    #print(zlevels[deep])
    td = (int(zlevels[deep]))
    #print(td)
    # AX 1 
    zmax, ylocn = 41, 647
    section_slice = np.arange(110, 230)
    #print('shape of plotdat')
    #print(plotdat.shape)
    pdat = plotdat[t,:,ylocn,section_slice]
    #print('shape of pdat')
    #print(pdat.shape)
    pdat = np.ma.masked_values(plotdat[t,:,ylocn,section_slice],0)
    pdat = pdat.T
    cmap = tcmap
    ss = section_slice[:]
    #print(zlevels.shape)
    #print(ss.shape)
    x, y, = np.meshgrid(section_slice,zlevels)
    #print(x.shape)
    
    mesh = ax1.pcolormesh(
        x, y, pdat, cmap=tcmap, vmin=v_min, vmax=v_max) #,
   #     ,
   # )
    ax1.invert_yaxis()

    ax1.set_xlabel('x Index', fontsize = 15)
    ax1.set_ylabel('depth (m)', fontsize = 15)
    # Axes limits and grid
    ax1.set_xlim(section_slice[1], section_slice[-1])
    ax1.set_ylim(zlevels[zmax - 2] + 20, 0)
    ax1.set_facecolor(land_colour)
    #ax1.grid()
    ax1.set_title('NORTH Crossection', fontsize = 20, color = 'xkcd:bright blue')

    # AX 2
    zmax, ylocn = 41, 505
    section_slice = np.arange(190, 310)
    #pdat = plotdat[t,:,ylocn,section_slice]
    pdat = np.ma.masked_values(plotdat[t,:,ylocn,section_slice],0)
    pdat = pdat.T
    cmap = tcmap
    mesh = ax2.pcolormesh(
        section_slice[:], zlevels[:zmax], pdat,
        cmap=cmap, vmin=v_min, vmax=v_max,
    )
    ax2.invert_yaxis()
    ax2.set_xlabel('x Index', fontsize = 18)
    # Axes limits and grid
    ax2.set_xlim(section_slice[1], section_slice[-1])
    ax2.set_ylim(zlevels[zmax - 2] + 20, 0)
    ax2.set_facecolor(land_colour)
    #ax2.grid()
    ax2.set_yticks([])
    ax2.set_title('NORTH - CENTRAL Crossection', fontsize = 20, color = 'xkcd:royal blue')
    
    # AX 3
    
    zmax, ylocn = 41, 423
    section_slice = np.arange(200, 320)
    #pdat = plotdat[t,:,ylocn,section_slice]
    pdat = np.ma.masked_values(plotdat[t,:,ylocn,section_slice],0)
    pdat = pdat.T
    cmap = tcmap

    mesh = ax3.pcolormesh(
        section_slice[:], zlevels[:zmax], pdat,
        cmap=cmap, vmin=v_min, vmax=v_max,
    )
    ax3.invert_yaxis()
    ax3.set_xlabel('x Index', fontsize = 18)
    # Axes limits and grid
    ax3.set_xlim(section_slice[1], section_slice[-1])
    ax3.set_ylim(zlevels[zmax - 2] + 20, 0)
    ax3.set_facecolor(land_colour)
    #ax3.grid()
    ax3.set_yticks([])
    ax3.set_title('SOUTH Crossection', fontsize = 20, color = 'xkcd:red')
    
    
    ### ax4 is the map
    
    #x_slice = np.arange(0, 398)
    #y_slice = np.arange(0, 898)

    surf_dat =  np.ma.masked_values(plotdat[t, deep, :, :], 0)

    viz_tools.set_aspect(ax4)
    
    line_s = np.arange(200, 320)
    ax4.plot(
        line_s, 423*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='red',
        label='Section Line',
    )
    i = 263
    j = 423
    ax4.plot(i,j,marker='o',color='xkcd:red',markersize=12, markeredgecolor = 'w', markeredgewidth = 2.0)
    ax4.text(i-150,j+2,'SOUTH',fontsize = 13)
    
    line_s = np.arange(190, 310)
    ax4.plot(
        line_s, 505*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='xkcd:royal blue',
        label='Section Line',
    )
    i = 246
    j = 505
    ax4.plot(i,j,marker='o',color='xkcd:royal blue',markersize=12, markeredgecolor = 'w', markeredgewidth = 2.0)
    ax4.text(i-250,j+2,'N.-CENTRAL',fontsize = 13)
    
    
    
    line_s = np.arange(110, 230)
    ax4.plot(
        line_s, 647*np.ones_like(line_s),
        linestyle='solid', linewidth=1.5, color='xkcd:bright blue',
        label='Section Line',
    )
    i = 168
    j = 647
    ax4.plot(i,j,marker='o',color='xkcd:bright blue',markersize=12, markeredgecolor = 'w', markeredgewidth = 2.0)
    print('text?')
    ax4.text(i-150,j+2,'NORTH',fontsize = 13)
    

    
    cmap.set_bad(land_colour)
    mesh = ax4.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)
    ax4.set_xlabel('')
    ax4.set_ylabel('')
    #ax4.set_xticks([])
    #ax4.set_yticks([])

    ax4.grid()
    cbar = fig.colorbar(mesh, ax=ax4)
    cbar.set_label(cblabel, fontsize = 15)
    cbar.set_ticks([])
    if td < 1:
        ax4.set_title('Surface ' + surtit, fontsize = 20, color = 'black')
    else:
        ax4.set_title(surtit +'at depth ' + str(td) +' m', fontsize = 20, color = 'black')
    #####


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
    
    tit = surtit + tit2
    plt.suptitle(tit, fontsize = 20)
    total_fig = dirstr+figtit+ lsi + '.png'
    fig.savefig(total_fig)
    #plt.show()
    plt.close(fig)


    
def MRAP(resdir,dirstr,segment,pdat,deep):

    day_indexer = [0,15,15,14,16,15,15,15]   
    indexer = (sum(day_indexer[0:segment]))
    print('index')
    print(indexer)
    
    trac, ptrac, vtrac, ratrac, timesteps  = dat_retrieve_DAY(resdir)
    print('number of frames')
    print(timesteps.shape)
    
    dic = trac.variables['dissolved_inorganic_carbon'][:]
    oxy = trac.variables['dissolved_oxygen'][:]
    yvel = vtrac.variables['vomecrty'][:]
    sn = ptrac.variables['vosaline'][:]
    


    if pdat == 'oxy':
        plotdat = oxy
        tcmap = cm.cm.speed
        cblabel = ''
        surtit = 'BWT ' 
        v_min = 7
        v_max = 15



    w = (timesteps.shape)

    timer = w[0]
    print(timer)
    for t in range(0,timer):
        
        figtit = pdat + '_zl_' + str(deep) + 'range'+ str(v_min) + str(v_max) + '_3s'
        print(figtit)
        a3stn_plot(plotdat, t, v_min, v_max, tcmap, dirstr, figtit, indexer, surtit, deep, cblabel)
        
        figtit = pdat + '_zl_' + str(deep) + 'range'+ str(v_min) + str(v_max) + '_sp'
        print(figtit)
        
        surf_plot(plotdat, t, v_min, v_max, tcmap, dirstr, figtit, indexer,surtit, deep, cblabel) 
        
