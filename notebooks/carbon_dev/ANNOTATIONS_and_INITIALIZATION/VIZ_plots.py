def two_panel_plot(surfdat_1,surfdat_2,tit1,tit2,t_cmap,xsize,ysize,v_min1,v_max1,v_min2,v_max2,cl1,cl2,bigtit):
    "TESTED"
    from salishsea_tools import (teos_tools, tidetools, viz_tools)
    import matplotlib.pyplot as plt
    import cmocean as cm
    import numpy as np
    fig, axs = plt.subplots(1, 2, figsize=(xsize, ysize), sharey=True)
    cmap = t_cmap
            
    time_steps = (0,1)
    for ax, t in zip(axs, time_steps):

        if t == 0 :
            tplt = np.ma.masked_values(surfdat_1,0)
            ax.set_title(tit1)
            v_min = v_min1
            v_max = v_max1
            clabel = cl1
        if t == 1 :
            tplt = np.ma.masked_values(surfdat_2,0)
            ax.set_title(tit2)
            v_min = v_min2
            v_max = v_max2
            clabel = cl2

        viz_tools.set_aspect(ax)
        

        mesh = ax.pcolormesh(tplt, cmap=t_cmap, vmin=v_min, vmax=v_max)
        
        cbar = fig.colorbar(mesh, ax=ax)
        
        cbar.set_label(clabel)
        ax.set_xlabel('x Index')

    axs[0].set_ylabel('y Index')
    
    
    cmap.set_bad('slategray')
    plt.suptitle(bigtit,fontsize=20)

def animate_surf(hrly_dat,tit,v_min,v_max,t_cmap):
    "TESTED"
    from salishsea_tools import (teos_tools, tidetools, viz_tools)
    import matplotlib.pyplot as plt
    import cmocean as cm
    import numpy as np
    from matplotlib import animation, rc
    rc('animation', html='html5')
    fig,ax = plt.subplots(figsize=(11,11))
    
    def init():
        tplt = hrly_dat[0,:,:]
        tplt = np.ma.masked_values(tplt, 0)  
        
        s = ax.pcolormesh(tplt, cmap=t_cmap, vmin = v_min, vmax = v_max)
        t_cmap.set_bad('slategray')
        ax.set_title(tit + ', hour = 0') 
        fig.colorbar(s, ax=ax)

 
    def make_plot(i):
        ax.clear()
        viz_tools.set_aspect(ax)
        cmap = t_cmap
        tplt = hrly_dat[i,:,:]
        tplt = np.ma.masked_values(tplt, 0)
        mesh = ax.pcolormesh(tplt, cmap=t_cmap, vmin = v_min, vmax = v_max)
        t_cmap.set_bad('slategray')
        ax.set_title(tit + ', hour = ' + str(i+1)) 
        #ax.colorbar(mesh)
        #fig.colorbar(mesh, ax=ax)
        return ax, mesh
    
    interval = 0.3
    ani40 = animation.FuncAnimation(fig,make_plot,frames=24,interval=interval*1e+3, init_func = init, repeat=True)

    return ani40



def nice_CN_plot(tit1, tit2, plotdat, physdat, t, v_min, v_max, tcmap, clabel):
    "NOT TESTED"
    import netCDF4 as nc
    import matplotlib.pyplot as plt
    import datetime
    import os
    import numpy as np
    import cmocean as cm
    from salishsea_tools import visualisations as vis
    from salishsea_tools import (teos_tools, tidetools, viz_tools)
    import cmocean
    import glob
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(16, 8 ))
    land_colour = 'whitesmoke'
    

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
    cbar.set_label(clabel)

    # Axes labels and title
    axl.set_xlabel('x Index')
    axl.set_ylabel('depth (m)')
    axl.set_title(tit1)

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
    

    #cmap.set_bad('whitesmoke')
    mesh = axr.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)

    # Axes labels and title
    axr.set_xlabel('')
    axr.set_ylabel('')
    axr.set_xticks([])
    axr.set_yticks([])
    axr.set_title(tit2)
    legend = axr.legend(loc='best', fancybox=True, framealpha=0.25)
    axr.grid()
    
    
def ONC_animator(plotdat,tstart,tend,v_min,v_max,stepsize,tit1,tit2,figtit,dirstr,tcmap,clabel,indexer):
    import numpy as np
    import netCDF4 as nc
    import cmocean as cm
    import matplotlib.pyplot as plt
    from salishsea_tools import (teos_tools, tidetools, viz_tools)
    
    
    for i in range(tstart,tend):
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(16, 8 ))
        land_colour = 'whitesmoke'

        daypath = '/results/SalishSea/nowcast-blue/01dec17/SalishSea_1d_20171201_20171201_grid_T.nc' 
        t_d = nc.Dataset(daypath)

        zlevels = t_d.variables['deptht']

        # Define the component slice to plot
        zmax, ylocn = 41, 424
        section_slice = np.arange(208, 293)

        pdat = np.ma.masked_values(plotdat[i,:,424,section_slice],0)

        cmap = tcmap
        mesh = axl.pcolormesh(
            section_slice[:], zlevels[:zmax], pdat,
            cmap=cmap, vmin=v_min, vmax=v_max,
        )
        axl.invert_yaxis()
        cbar = fig.colorbar(mesh, ax=axl)
        cbar.set_label(clabel)

        # Axes labels and title
        axl.set_xlabel('x Index')
        axl.set_ylabel('depth (m)')
        axl.set_title(tit1)

        # Axes limits and grid
        axl.set_xlim(section_slice[1], section_slice[-1])
        axl.set_ylim(zlevels[zmax - 2] + 10, 0)
        axl.set_facecolor(land_colour)
        axl.grid()

        # Define surface current magnitude slice
        x_slice = np.arange(0, 398)
        y_slice = np.arange(0, 898)
        line_s = np.arange(0,398)
        
        surf_dat =  np.ma.masked_values(plotdat[i, 0, y_slice, x_slice], 0)
        viz_tools.set_aspect(axr)
        axr.plot(
            line_s, 424*np.ones_like(line_s),
            linestyle='solid', linewidth=3, color='black',
            label='Section Line',
        )

        cmap.set_bad(land_colour)
        mesh = axr.pcolormesh(surf_dat, cmap=tcmap, vmin = v_min, vmax = v_max)

        # Axes labels and title
        axr.set_xlabel('')
        axr.set_ylabel('')
        axr.set_xticks([])
        axr.set_yticks([])
        axr.set_title(tit2)
        legend = axr.legend(loc='best', fancybox=True, framealpha=0.25)
        axr.grid()
        
        t_i = indexer + i
        si = str(t_i)
        #print(len(si))
        if len(si) == 1:
            lsi = '00' + si
        if len(si) == 2:
            lsi = '0' + si
        if len(si) == 3:
            lsi = si
        print(lsi)
        fig.suptitle(figtit + str(i))
        fig.savefig(dirstr+figtit+str(lsi))
        plt.close(fig)


def thalweg_animator(plotdat,tstart,tend,vmin,vmax,stepsize,tit,tit2,figtit,dirstr,indexer,t_cmap,clabel):
    "TESTED"
    from salishsea_tools import visualisations as vis
    import cmocean as cm
    import matplotlib.pyplot as plt
    import netCDF4 as nc
    import numpy as np 
    
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
            lsi = '00' + si
        if len(si) == 2:
            lsi = '0' + si
        if len(si) == 3:
            lsi = si
        print(lsi)
        fig.suptitle(tit + str(lsi) + tit2, fontsize = 18)
        fig.savefig(dirstr+figtit+str(lsi))
        plt.close(fig)
        
        
def phpco2om_plot(surfdat_1,surfdat_2,surfdat_3,tit1,tit2,tit3,t_cmap1,t_cmap2,t_cmap3,xsize,ysize,v_min1,v_max1,v_min2,v_max2,v_min3,v_max3,cl1,cl2,cl3,bigtit,dirstr,ind):
    "TESTED"
    from salishsea_tools import (teos_tools, tidetools, viz_tools)
    import matplotlib.pyplot as plt
    import cmocean as cm
    import numpy as np
    fig, axs = plt.subplots(1, 3, figsize=(xsize, ysize), sharey=True)
    
            
    time_steps = (0,1,2)
    for ax, t in zip(axs, time_steps):

        if t == 0 :
            tplt = np.ma.masked_values(surfdat_1,0)
            ax.set_title(tit1)
            v_min = v_min1
            v_max = v_max1
            t_cmap = t_cmap1
            clabel = cl1
        if t == 1 :
            tplt = np.ma.masked_values(surfdat_2,0)
            ax.set_title(tit2)
            v_min = v_min2
            v_max = v_max2
            t_cmap = t_cmap2
            clabel = cl2

        if t == 2 :
            tplt = np.ma.masked_values(surfdat_3,0)
            ax.set_title(tit3)
            v_min = v_min3
            v_max = v_max3
            t_cmap = t_cmap3
            clabel = cl3
        viz_tools.set_aspect(ax)
        

        mesh = ax.pcolormesh(tplt, cmap=t_cmap, vmin=v_min, vmax=v_max)
        
        cbar = fig.colorbar(mesh, ax=ax)
        
        cbar.set_label(clabel)
        ax.set_xlabel('x Index')

    axs[0].set_ylabel('y Index')
    
    
    t_cmap.set_bad('slategray')
    plt.suptitle(bigtit,fontsize=20)
    figtit = 'CCpar_day_' + str(ind)
    total_fig = dirstr+figtit + '.png'
    fig.savefig(total_fig)
    plt.close(fig)
    
    
    
def one_panel_plot(surfdat_1,tit1,t_cmap,xsize,ysize,v_min1,v_max1,cl1,bigtit):
    "TESTED"
    from salishsea_tools import (teos_tools, tidetools, viz_tools)
    import matplotlib.pyplot as plt
    import cmocean as cm
    import numpy as np
    fig, ax = plt.subplots(1, 1, figsize=(xsize, ysize), sharey=True)
    cmap = t_cmap
            

    tplt = np.ma.masked_values(surfdat_1,0)
    ax.set_title(tit1)
    v_min = v_min1
    v_max = v_max1
    clabel = cl1

    viz_tools.set_aspect(ax)


    mesh = ax.pcolormesh(tplt, cmap=t_cmap, vmin=v_min, vmax=v_max)

    cbar = fig.colorbar(mesh, ax=ax)

    cbar.set_label(clabel)
    ax.set_xlabel('x Index')

    ax.set_ylabel('y Index')
    
    
    cmap.set_bad('slategray')
    plt.suptitle(bigtit,fontsize=20)
        