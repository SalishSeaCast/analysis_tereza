import arrow
import netCDF4 as nc
import CCCma_profiles as pf
import CCCma_stations as cs
import CCCma_ncmaker as nm
import numpy as np
import gsw
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import reload
import cmocean as cm

from salishsea_tools import (
#    nc_tools,
    viz_tools
#    geo_tools,
#    tidetools
)

def parameterspace4(pars_profs,stns,rdir,ddmmmyy,humandate, dss_sig):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=(19, 12) , nrows=2, ncols=2)

    fs = 18
    cols = []
    xs = []
    ys = []
    stn_in = []
    for s in stns:
        col = stns[s]['color']
        x = stns[s]['x']
        y = stns[s]['y']
        stn = stns[s]['code']
        cols.append(col)
        xs.append(x)
        ys.append(y)
        stn_in.append(stn)
    nos = len(cols)
    
    for i in range(0,nos):
        sal_prof = pars_profs['sal']
        DIC_prof = pars_profs['DIC']
        ax1.plot(sal_prof[i,:],DIC_prof[i,:],c=cols[i],marker='o',
                           markersize=15, linestyle = '')
        
        ax1.set_title("DIC vs Salinity",fontsize=fs)
    #TA
    for i in range(0,nos):
        TA_prof = pars_profs['TA']
        ax2.plot(sal_prof[i,:],TA_prof[i,:],c=cols[i],marker='o',
                           markersize=15, linestyle = '')
        ax2.set_title("TA vs Salinity",fontsize=fs)
        ax2.legend(stn_in , fontsize = fs, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #OmegaA
    for i in range(0,nos):
        OmA_prof = pars_profs['OmA']
        ax3.plot(sal_prof[i,:],OmA_prof[i,:],c=cols[i],marker='o',
                           markersize=15, linestyle = '')
        ax3.set_title("OmegaA vs Salinity",fontsize=fs)
        
    #OmegaA
    for i in range(0,nos):
        pH_prof = pars_profs['pH']
        ax4.plot(sal_prof[i,:],pH_prof[i,:],c=cols[i],marker='o',
                           markersize=15, linestyle = '')
        ax4.set_title("pH vs Salinity",fontsize=fs)

    ax1.set_xlabel('Salinity (g/kg)', fontsize = fs)
    ax2.set_xlabel('Salinity (g/kg)', fontsize = fs)
    ax3.set_xlabel('Salinity (g/kg)', fontsize = fs)
    ax4.set_xlabel('Salinity (g/kg)', fontsize = fs)
    
    ax1.tick_params(axis='y', labelsize=fs )
    ax2.tick_params(axis='y', labelsize=fs )
    ax3.tick_params(axis='y', labelsize=fs )
    ax4.tick_params(axis='y', labelsize=fs )

    ax1.tick_params(axis='x', labelsize=fs )
    ax2.tick_params(axis='x', labelsize=fs )
    ax3.tick_params(axis='x', labelsize=fs )
    ax4.tick_params(axis='x', labelsize=fs )
    
    ax1.set_xlim(15,35)
    ax2.set_xlim(15,35)
    ax3.set_xlim(15,35)
    ax4.set_xlim(15,35)
    
    ax1.set_ylim(1800,2300)
    ax2.set_ylim(1800,2300)
    ax3.set_ylim(0,2)
    ax4.set_ylim(6.8,8.2)
    
    st = 'Salish Sea Carbonate Chemistry Parameterspace Diagram, ' + humandate
    plt.suptitle(st,fontsize=fs)
    fname = rdir + f'{ddmmmyy}_4pspace_' + dss_sig +'.png'

    fig.savefig(fname)
    plt.close()