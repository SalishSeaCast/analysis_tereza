#one_panel_plot
#profiles
#point_value
#profile_plotter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.style.use('seaborn-whitegrid')
import netCDF4 as nc
import cmocean as cm
import numpy as np
from salishsea_tools import (
    viz_tools,
)
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
#import CCCma
import CCCma_stations as cs
#from matplotlib import reload
import arrow
import gsw


def one_panel_plot(surfdat_1, stns, tit1,t_cmap,xsize,ysize,v_min1,v_max1,cl1,bigtit):
    "TESTED"

    fig, ax = plt.subplots(1, 1, figsize=(xsize, ysize), sharey=True)
    cmap = t_cmap
    tplt = np.ma.masked_values(surfdat_1,0)
    ax.set_title(tit1,fontsize = 20 )
    v_min = v_min1
    v_max = v_max1
    clabel = cl1

    viz_tools.set_aspect(ax)
    mesh = ax.pcolormesh(tplt, cmap=t_cmap, vmin=v_min, vmax=v_max)
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(clabel, fontsize = 20 )
    cbar.ax.tick_params(labelsize=20)
 
    ax.set_xlabel('x Index', fontsize = 20 )
    ax.set_ylabel('y Index', fontsize = 20 )
    
    for key in stns:

        tx = stns[key]['x']
        ty = stns[key]['y']
        col = stns[key]['color']
        code = stns[key]['code']
        pat = patches.Rectangle((tx,ty),20,20,linewidth=2,edgecolor=col,facecolor='none')
        ax.add_patch(pat)
        ax.text(tx+22,ty+3,code, weight = 'bold', fontsize = 20)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    cmap.set_bad('aliceblue')
    plt.suptitle(bigtit,fontsize=20)
    
def profiles(carp, grid, stns):
             
    '''
    Take 2 daily datasets, carbon+ and grid, extract depth 
    profiles of sal,temp,DIC,TA,O2,pH,OmA and their standard deviations"
    '''
    w = carp
    wp = grid

    
    sal = wp.variables['vosaline'][0,:,:,:]
    temp = wp.variables['votemper'][0,:,:,:]
    DIC = w.variables['dissolved_inorganic_carbon'][0,:,:,:]
    TA = w.variables['total_alkalinity'][0,:,:,:]
    O2 = w.variables['dissolved_oxygen'][0,:,:,:]
    prof_depth = w.variables['deptht'][:]

    
    depth_broadcast = np.zeros([40,20,20])
    for i in range(0,40):
        depth_broadcast[i,:,:] = prof_depth[i]

    nos = len(stns)
    stn_list = []
    
    sal_prof = np.zeros([nos,40])
    temp_prof = np.zeros([nos,40])
    DIC_prof = np.zeros([nos,40])
    TA_prof = np.zeros([nos,40])
    pH_prof = np.zeros([nos,40])
    OmA_prof = np.zeros([nos,40])
    O2_prof = np.zeros([nos,40])

    sal_profSD = np.zeros([nos,40])
    temp_profSD = np.zeros([nos,40])
    DIC_profSD = np.zeros([nos,40])
    TA_profSD = np.zeros([nos,40])
    pH_profSD = np.zeros([nos,40])
    OmA_profSD = np.zeros([nos,40])
    O2_profSD = np.zeros([nos,40])
    
    b = 0
    for s in stns:
        stn_list.append(s)
        stn = stns[s]
        #print('Calculating profiles for ' + stns[s]['fullname'])
        ty = stns[s]['y']
        tx = stns[s]['x']

        ts = sal[:,ty:ty+20,tx:tx+20]
        tte = temp[:,ty:ty+20,tx:tx+20]
        td = DIC[:,ty:ty+20,tx:tx+20]
        tta = TA[:,ty:ty+20,tx:tx+20]
        to2 = O2[:,ty:ty+20,tx:tx+20]

        tsr =  ts.reshape(40,400)
        tter =  tte.reshape(40,400)
        tdr =  td.reshape(40,400)
        ttar =  tta.reshape(40,400)
        tdepth = depth_broadcast.reshape(40,400)
        to2r = to2.reshape(40,400)

        tsr[tsr == 0] = np.ma.masked
        tter[tter == 0] = np.ma.masked
        tdr[tdr == 0] = np.ma.masked
        ttar[ttar == 0] = np.ma.masked

        sal_prof[b,:] = np.ma.MaskedArray.nanmean(tsr, axis = 1)
        sal_profSD[b,:] = np.ma.MaskedArray.nanstd(tsr, axis = 1)
        
        temp_prof[b,:] = np.ma.MaskedArray.nanmean(tter, axis = 1)
        temp_profSD[b,:] = np.ma.MaskedArray.nanstd(tter, axis = 1)

        DIC_prof[b,:] = np.ma.MaskedArray.nanmean(tdr, axis = 1)
        DIC_profSD[b,:] = np.ma.MaskedArray.nanstd(tdr, axis = 1)

        TA_prof[b,:] = np.ma.MaskedArray.nanmean(ttar, axis = 1)
        TA_profSD[b,:] = np.ma.MaskedArray.nanstd(ttar, axis = 1)

        O2_prof[b,:] = np.ma.MaskedArray.nanmean(to2r, axis = 1)
        O2_profSD[b,:] = np.ma.MaskedArray.nanstd(to2r, axis = 1)

        tsra = np.ravel(ts)
        ttera = np.ravel(tte)
        tdra = np.ravel(td) * 1e-3
        ttara =  np.ravel(tta) * 1e-3
        tdepthra = np.ravel(tdepth)
        tpressure = tdepthra 
        tpressure[:] =1
        tzero = tdepthra * 0 
        
        tsra_psu = tsra*35/35.16504
        ttera_is = gsw.t_from_CT(tsra,ttera,tdepthra)

        response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                           sil=tzero, phos=tzero, patm=tpressure, depth=tdepthra, lat=tzero, 
                            optcon='mol/m3', optt='Tinsitu', optp='m',
                            optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
        pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

        pHra = pH.reshape(40,400)
        OmegaAra = OmegaA.reshape(40,400)

        pHra[pHra == 1.00000000e+20] = np.nan
        OmegaAra[OmegaAra == 1.00000000e+20] = np.nan
        pHra = np.ma.masked_invalid(pHra)
        OmegaAra = np.ma.masked_invalid(OmegaAra)

        pH_prof[b,:] = np.ma.MaskedArray.nanmean(pHra, axis = 1)
        pH_profSD[b,:] = np.ma.MaskedArray.nanstd(pHra, axis = 1)

        OmA_prof[b,:] = np.ma.MaskedArray.nanmean(OmegaAra, axis = 1)
        OmA_profSD[b,:] = np.ma.MaskedArray.nanstd(OmegaAra, axis = 1)
        
        #depth = w.variables['deptht'][:]
        sal_prof[sal_prof == 0 ] = np.nan
        temp_prof[temp_prof == 0 ] = np.nan
        DIC_prof[DIC_prof == 0 ] = np.nan
        TA_prof[TA_prof == 0 ] = np.nan
        pH_prof[pH_prof == 0 ] = np.nan
        OmA_prof[OmA_prof == 0 ] = np.nan
        
        b = b+1
        
    pars_profs = {'sal': sal_prof, 'sal_SD': sal_profSD,\
         'temp': temp_prof, 'temp_SD': temp_profSD,\
         'DIC': DIC_prof, 'DIC_SD' : DIC_profSD,\
         'TA': TA_prof, 'TA_SD' : TA_profSD,\
         'OmA': OmA_prof, 'OmA_SD': OmA_profSD,\
         'pH': pH_prof, 'pH_SD': pH_profSD,\
         'O2' : O2_prof, 'O2_SD': O2_profSD}
    return pars_profs, stn_list, prof_depth

def point_value(carp, grid, stns):
    '''
    
    '''   
    w = carp
    wp = grid
    
    sal = wp.variables['vosaline'][0,:,:,:]
    temp = wp.variables['votemper'][0,:,:,:]
    DIC = w.variables['dissolved_inorganic_carbon'][0,:,:,:]
    TA = w.variables['total_alkalinity'][0,:,:,:]
    O2 = w.variables['dissolved_oxygen'][0,:,:,:]
    depth = w.variables['deptht'][:]
    
    depth_broadcast = np.zeros([40,20,20])
    for i in range(0,40):
        depth_broadcast[i,:,:] = depth[i]
    
    depth_inds = [0,21,26]
    dum = [0.0,0.0,0.0]
    pt_depths = np.zeros_like(dum)
    
    for i in range(0,len(pt_depths)):
        di = depth_inds
        pt_depths[i] = depth[depth_inds[i]]
    
    nos = len(stns)
    stn_list = []

    sal_pt = np.zeros([nos,len(depth_inds)])
    temp_pt = np.zeros([nos,len(depth_inds)])
    DIC_pt = np.zeros([nos,len(depth_inds)])
    TA_pt = np.zeros([nos,len(depth_inds)])
    pH_pt = np.zeros([nos,len(depth_inds)])
    OmA_pt = np.zeros([nos,len(depth_inds)])
    O2_pt = np.zeros([nos,len(depth_inds)])

    sal_ptSD = np.zeros([nos,len(depth_inds)])
    temp_ptSD = np.zeros([nos,len(depth_inds)])
    DIC_ptSD = np.zeros([nos,len(depth_inds)])
    TA_ptSD = np.zeros([nos,len(depth_inds)])
    pH_ptSD = np.zeros([nos,len(depth_inds)])
    OmA_ptSD = np.zeros([nos,len(depth_inds)])
    O2_ptSD = np.zeros([nos,len(depth_inds)])
    
    b = 0
    for s in stns:
        stn_list.append(s)
        
        
        #print('Calculating point values for ' + stns[s]['fullname'])
        ty = stns[s]['y']
        tx = stns[s]['x']
        
        for d in range(0,len(depth_inds)):
            dp = depth[depth_inds[d]]
            ts = sal[d,ty:ty+20,tx:tx+20]
            tte = temp[d,ty:ty+20,tx:tx+20]
            td = DIC[d,ty:ty+20,tx:tx+20]
            tta = TA[d,ty:ty+20,tx:tx+20]
            to2 = O2[d,ty:ty+20,tx:tx+20]

            tsr =  ts.reshape(1,400)
            tter =  tte.reshape(1,400)
            tdr =  td.reshape(1,400)
            ttar =  tta.reshape(1,400)

            tdepth = np.zeros_like(ttar)
     

            tdepth[:] = dp
            to2r = to2.reshape(1,400)

            tsr[tsr == 0] = np.ma.masked
            tter[tter == 0] = np.ma.masked
            tdr[tdr == 0] = np.ma.masked
            ttar[ttar == 0] = np.ma.masked

            sal_pt[b,d] = np.ma.MaskedArray.nanmean(tsr, axis = 1)
            sal_ptSD[b,d] = np.ma.MaskedArray.nanstd(tsr, axis = 1)

            temp_pt[b,:] = np.ma.MaskedArray.nanmean(tter, axis = 1)
            temp_ptSD[b,d] = np.ma.MaskedArray.nanstd(tter, axis = 1)

            DIC_pt[b,d] = np.ma.MaskedArray.nanmean(tdr, axis = 1)
            DIC_ptSD[b,d] = np.ma.MaskedArray.nanstd(tdr, axis = 1)

            TA_pt[b,d] = np.ma.MaskedArray.nanmean(ttar, axis = 1)
            TA_ptSD[b,d] = np.ma.MaskedArray.nanstd(ttar, axis = 1)

            O2_pt[b,d] = np.ma.MaskedArray.nanmean(to2r, axis = 1)
            O2_ptSD[b,d] = np.ma.MaskedArray.nanstd(to2r, axis = 1)

            tsra = np.ravel(ts)
            ttera = np.ravel(tte)
            tdra = np.ravel(td) * 1e-3
            ttara =  np.ravel(tta) * 1e-3
            tdepthra = np.ravel(tdepth)
            tpressure = np.zeros_like(tdepthra)
            tpressure[:] =1
            tzero = tdepthra * 0 
            
            tsra_psu = tsra*35/35.16504
            ttera_is = gsw.t_from_CT(tsra,ttera,tdepthra)

            response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                               sil=tzero, phos=tzero, patm=tpressure, depth=tdepthra, lat=tzero, 
                                optcon='mol/m3', optt='Tinsitu', optp='m',
                                optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
            pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

            pHra = pH.reshape(1,400)
            OmegaAra = OmegaA.reshape(1,400)

            pHra[pHra == 1.00000000e+20] = np.nan
            OmegaAra[OmegaAra == 1.00000000e+20] = np.nan
            pHra = np.ma.masked_invalid(pHra)
            OmegaAra = np.ma.masked_invalid(OmegaAra)
        
            pH_pt[b,d] = np.ma.MaskedArray.nanmean(pHra, axis = 1)
            pH_pt[b,d] = np.ma.MaskedArray.nanstd(pHra, axis = 1)

            OmA_pt[b,d] = np.ma.MaskedArray.nanmean(OmegaAra, axis = 1)
            OmA_pt[b,d] = np.ma.MaskedArray.nanstd(OmegaAra, axis = 1)
            
            sal_pt[sal_pt == 0 ] = np.nan
            temp_pt[temp_pt == 0 ] = np.nan
            DIC_pt[DIC_pt == 0 ] = np.nan
            TA_pt[TA_pt == 0 ] = np.nan
            pH_pt[pH_pt == 0 ] = np.nan
            OmA_pt[OmA_pt == 0 ] = np.nan
            
        wb = b+1

    pars_pts = {'sal': sal_pt, 'sal_SD': sal_ptSD,\
     'temp': temp_pt, 'temp_SD': temp_ptSD,\
     'DIC': DIC_pt, 'DIC_SD' : DIC_ptSD,\
     'TA': TA_pt, 'TA_SD' : TA_ptSD,\
     'OmA': OmA_pt, 'OmA_SD': OmA_ptSD,\
     'pH': pH_pt, 'pH_SD': pH_ptSD,\
     'O2' : O2_pt, 'O2_SD': O2_ptSD}
    
    return pars_pts, pt_depths, stn_list    
    
def profile_plotter(pars_profs,depth,stns,date,ddmmmyy, rdir, dss_sig):
    fa = 1

    fig, ((ax1, ax2, ax3, ax4, ax5, ax6), (ax7, ax8, ax9, ax10, ax11, ax12)) = \
    plt.subplots(figsize=(14*fa, 12.0*fa) , nrows=2, ncols=6)

    nos = len(stns)

    cols = []
    codes = []
    for s in stns:

        col = stns[s]['color']
        code = stns[s]['code']
        cols.append(col)
        codes.append(code)

    for i in range(0,nos):

        ## salinity
        sal_prof = pars_profs['sal']
        sal_profSD = pars_profs['sal_SD']
        pattern = sal_prof[i,:]
        st_dev = sal_profSD[i,:]
        x1 = pattern-st_dev
        x2 = pattern + st_dev
        ax7.plot(sal_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax7.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax7.set_ylabel('Depth (m)', fontsize = 20)
        ax7.set_ylim([0,350])
        ax7.set_ylim(ax7.get_ylim()[::-1])
        ax1.plot(sal_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i])
        ax1.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax1.set_ylabel('Depth (m)', fontsize = 20)
        ax1.set_ylim([0,20])
        ax1.set_ylim(ax1.get_ylim()[::-1])

         ## temp
        temp_prof = pars_profs['temp']
        temp_profSD = pars_profs['temp_SD']
        pattern = temp_prof[i,:]
        st_dev = temp_profSD[i,:]
        x1 = pattern-st_dev
        x2 = pattern + st_dev
        ax8.plot(temp_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax8.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax8.set_ylim([0,350])
        ax8.set_ylim(ax8.get_ylim()[::-1])
        ax2.plot(temp_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax2.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax2.set_ylim([0,20])
        ax2.set_ylim(ax2.get_ylim()[::-1])

         ## DIC
        DIC_prof = pars_profs['DIC']
        DIC_profSD = pars_profs['DIC_SD']
        pattern = DIC_prof[i,:]
        st_dev = DIC_profSD[i,:]
        x1 = pattern-st_dev
        x2 = pattern + st_dev
        ax9.plot(DIC_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax9.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax9.set_ylim([0,350])
        ax9.set_ylim(ax9.get_ylim()[::-1])
        ax3.plot(DIC_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax3.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax3.set_ylim([0,20])
        ax3.set_ylim(ax3.get_ylim()[::-1])

    #     ## TA
        TA_prof = pars_profs['TA']
        TA_profSD = pars_profs['TA_SD']
        pattern = TA_prof[i,:]
        st_dev = TA_profSD[i,:]
        x1 = pattern-st_dev
        x2 = pattern + st_dev
        ax10.plot(TA_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax10.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax10.set_ylim([0,350])
        ax10.set_ylim(ax10.get_ylim()[::-1])
        ax4.plot(TA_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax4.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax4.set_ylim([0,20])
        ax4.set_ylim(ax4.get_ylim()[::-1])

    #     ## ph
        pH_prof = pars_profs['pH']
        pH_profSD = pars_profs['pH_SD']
        pattern = pH_prof[i,:]
        st_dev = pH_profSD[i,:]
        x1 = pattern-st_dev
        x2 = pattern + st_dev
        ax11.plot(pH_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax11.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax11.set_ylim([0,350])
        ax11.set_ylim(ax11.get_ylim()[::-1])
        ax5.plot(pH_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax5.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax5.set_ylim([0,20])
        ax5.set_ylim(ax5.get_ylim()[::-1])

    #     ## oma
        OmA_prof = pars_profs['OmA']
        OmA_profSD = pars_profs['OmA_SD']
        pattern = OmA_prof[i,:]
        st_dev = OmA_profSD[i,:]
        x1 = pattern-st_dev
        x2 = pattern + st_dev
        ax12.plot(OmA_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax12.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax12.set_ylim([0,350])
        ax12.set_ylim(ax12.get_ylim()[::-1])
        ax6.plot(OmA_prof[i,:],depth, linestyle='-', linewidth=3, color = cols[i]) 
        ax6.fill_betweenx(depth, x1, x2, facecolor=cols[i], alpha = 0.2)
        ax6.set_ylim([0,20])
        ax6.set_ylim(ax6.get_ylim()[::-1])

    #ax1.legend(codes , fontsize = 20)
    ax6.legend(codes , fontsize = 20, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax1.set_ylabel('Near-Surface Values: Depth (m)', fontsize = 25)
    ax7.set_ylabel('Full Profiles: Depth (m)', fontsize = 25)

    #limits for x axis
    ax1.set_xlim([15,32])
    ax7.set_xlim([15,32])
    ax2.set_xlim([0,20])
    ax8.set_xlim([0,20])
    ax3.set_xlim([700,2300])
    ax9.set_xlim([700,2300])
    ax4.set_xlim([700,2300])
    ax10.set_xlim([700,2300])
    ax5.set_xlim([7,9])
    ax11.set_xlim([7,9])
    ax6.set_xlim([0,2])
    ax12.set_xlim([0,2])
    
    fst = 16
    ax7.tick_params(axis='x', labelsize=fst )
    ax8.tick_params(axis='x', labelsize=fst )
    ax9.tick_params(axis='x', labelsize=fst)
    ax10.tick_params(axis='x', labelsize=fst)
    ax11.tick_params(axis='x', labelsize=fst)
    ax12.tick_params(axis='x', labelsize=fst)
    ax1.tick_params(axis='y', labelsize=fst)
    ax7.tick_params(axis='y', labelsize=fst)

    ax1.set_xticklabels('')
    ax2.set_xticklabels('')
    ax3.set_xticklabels('')
    ax4.set_xticklabels('')
    ax5.set_xticklabels('')
    ax6.set_xticklabels('')
    ax2.set_yticklabels('')
    ax3.set_yticklabels('')
    ax4.set_yticklabels('')
    ax5.set_yticklabels('')
    ax6.set_yticklabels('')
    ax8.set_yticklabels('')
    ax9.set_yticklabels('')
    ax10.set_yticklabels('')
    ax11.set_yticklabels('')
    ax12.set_yticklabels('')

    fsl = 20
    ax7.set_xlabel('Sal. (g/kg)', fontsize = fsl)
    ax8.set_xlabel('Temp. ( °C)', fontsize = fsl)
    ax9.set_xlabel('DIC μmol/kg', fontsize = fsl)
    ax10.set_xlabel('TA μmol/kg', fontsize = fsl)
    ax11.set_xlabel('pH', fontsize = fsl)
    ax12.set_xlabel('ΩA', fontsize = fsl)
    
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    st = 'Salish Sea Carbonate Chemistry Snapshot, ' + date
    fig.suptitle(st, fontsize = 30)
    
    fname = rdir  + f'{ddmmmyy}_prof_' + dss_sig +'.png'

    fig.savefig(fname)
    plt.close()


    
def ncmaker(stn_list, depths, pt_depths, pars_profs, pars_pts, ddmmmyy, rdir):
    #ncname = 'sample.nc'
    ncname = rdir + f'{ddmmmyy}_prof.nc'
    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('Profiles')
    g.createDimension('stn',len(stn_list))
    g.createDimension('z', len(depths))
    g.createDimension('time', None)
    #

    depth = g.createVariable('depth','f4','z')
    depth[:] = depths

    stn = g.createVariable('station','S4','stn')
    for i in range(0,len(stn_list)):
        stn[i] = stn_list[i]
        
    sal_prof = g.createVariable('sal_prof','f4',('stn','z'))
    sal_prof[:] = pars_profs['sal'][:]
    temp_prof = g.createVariable('temp_prof','f4',('stn','z'))
    temp_prof[:] = pars_profs['temp'][:]
    DIC_prof = g.createVariable('DIC_prof','f4',('stn','z'))
    DIC_prof[:] = pars_profs['DIC'][:]
    TA_prof = g.createVariable('TA_prof','f4',('stn','z'))
    TA_prof[:] = pars_profs['TA'][:]
    OmA_prof = g.createVariable('OmA_prof','f4',('stn','z'))
    OmA_prof[:] = pars_profs['OmA'][:]
    pH_prof = g.createVariable('pH_prof','f4',('stn','z'))
    pH_prof[:] = pars_profs['pH'][:]
    O2_prof = g.createVariable('O2_prof','f4',('stn','z'))
    O2_prof[:] = pars_profs['O2'][:]

    sal_profSD = g.createVariable('sal_profSD','f4',('stn','z'))
    sal_profSD[:] = pars_profs['sal_SD'][:]
    temp_profSD = g.createVariable('temp_profSD','f4',('stn','z'))
    temp_profSD[:] = pars_profs['temp_SD'][:]
    DIC_profSD = g.createVariable('DIC_profSD','f4',('stn','z'))
    DIC_profSD[:] = pars_profs['DIC_SD'][:]
    TA_profSD = g.createVariable('TA_profSD','f4',('stn','z'))
    TA_profSD[:] = pars_profs['TA_SD'][:]
    OmA_profSD = g.createVariable('OmA_profSD','f4',('stn','z'))
    OmA_profSD[:] = pars_profs['OmA_SD'][:]
    pH_profSD = g.createVariable('pH_profSD','f4',('stn','z'))
    pH_profSD[:] = pars_profs['pH_SD'][:]
    O2_profSD = g.createVariable('O2_profSD','f4',('stn','z'))
    O2_profSD[:] = pars_profs['O2_SD'][:]

    h = f.createGroup('Station_Averages')
    h.createDimension('stn',len(stn_list))
    h.createDimension('z_pt', len(pt_depths))
    h.createDimension('time', None)

    pt_depth = h.createVariable('depth','f4','z_pt')
    pt_depth[:] = pt_depths

    stn = h.createVariable('station','S4','stn')
    for i in range(0,len(stn_list)):
        stn[i] = stn_list[i]
        
    sal_pt = h.createVariable('sal_pt','f4',('stn','z_pt'))
    sal_pt[:] = pars_pts['sal'][:]
    temp_pt = h.createVariable('temp_pt','f4',('stn','z_pt'))
    temp_pt[:] = pars_pts['temp'][:]
    DIC_pt = h.createVariable('DIC_pt','f4',('stn','z_pt'))
    DIC_pt[:] = pars_pts['DIC'][:]
    TA_pt = h.createVariable('TA_pt','f4',('stn','z_pt'))
    TA_pt[:] = pars_pts['TA'][:]
    OmA_pt = h.createVariable('OmA_pt','f4',('stn','z_pt'))
    OmA_pt[:] = pars_pts['OmA'][:]
    pH_pt = h.createVariable('pH_pt','f4',('stn','z_pt'))
    pH_pt[:] = pars_pts['pH'][:]
    O2_pt = h.createVariable('O2_pt','f4',('stn','z_pt'))
    O2_pt[:] = pars_pts['O2'][:]

    sal_ptSD = h.createVariable('sal_ptSD','f4',('stn','z_pt'))
    sal_ptSD[:] = pars_pts['sal_SD'][:]
    temp_ptSD = h.createVariable('temp_ptSD','f4',('stn','z_pt'))
    temp_ptSD[:] = pars_pts['temp_SD'][:]
    DIC_ptSD = h.createVariable('DIC_ptSD','f4',('stn','z_pt'))
    DIC_ptSD[:] = pars_pts['DIC_SD'][:]
    TA_ptSD = h.createVariable('TA_ptSD','f4',('stn','z_pt'))
    TA_ptSD[:] = pars_pts['TA_SD'][:]
    OmA_ptSD = h.createVariable('OmA_ptSD','f4',('stn','z_pt'))
    OmA_ptSD[:] = pars_pts['OmA_SD'][:]
    pH_ptSD = h.createVariable('pH_ptSD','f4',('stn','z_pt'))
    pH_ptSD[:] = pars_pts['pH_SD'][:]
    O2_ptSD = h.createVariable('O2_ptSD','f4',('stn','z_pt'))
    O2_ptSD[:] = pars_pts['O2_SD'][:]

    f.close()
    
def surface_maps(carp,grid,stns, ddmmmyy, rdir,humandate, dss_sig):
    
    tsal = grid.variables['vosaline'][0,0,:,:]
    ttemp = grid.variables['votemper'][0,0,:,:]
    tdic = carp.variables['dissolved_inorganic_carbon'][0,0,:,:]
    tta = carp.variables['total_alkalinity'][0,0,:,:]

    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    tpressure[:] =1
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504
    ttera_is = gsw.t_from_CT(tsra,ttera,tzero)

    response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tzero, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    pHr = pH.reshape(898,398)
    OmA = OmegaA.reshape(898,398)
    
    surf_dat = [tsal, tdic, tta, ttemp, pHr, OmA]
    
    vmins = [25,1800,1800,5,7.5,0]
    vmaxs = [32,2200,2200,15,8.5,2]
    msk = [0,0,0,0,1e20,1e20]
    cl = ['salinity psu', 'DIC umol/kg', 'TA umol/kg', 'temp deg C', 'pH', 'Omega A']
    t_cmap = [cm.cm.haline, cm.cm.matter, cm.cm.matter, cm.cm.thermal, cm.cm.speed, cm.cm.curl]

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = \
    plt.subplots(figsize=(20, 27) , nrows=2, ncols=3)

    viz_tools.set_aspect(ax1)
    viz_tools.set_aspect(ax2)
    viz_tools.set_aspect(ax3)
    viz_tools.set_aspect(ax4)
    viz_tools.set_aspect(ax5)
    viz_tools.set_aspect(ax6)

    i = 0
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax1.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax1)
    cbar.set_label(cl[i], fontsize = 20 )

    i = 1
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax2.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax2)
    cbar.set_label(cl[i], fontsize = 20 )

    i = 2
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax3.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax3)
    cbar.set_label(cl[i], fontsize = 20 )

    i = 3
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax4.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax4)
    cbar.set_label(cl[i], fontsize = 20 )

    i = 4
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax5.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax5)
    cbar.set_label(cl[i], fontsize = 20 )

    i = 5
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax6.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax6)
    cbar.set_label(cl[i], fontsize = 20 )

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
        
    for w in range(0,len(stns)):        
        pat = patches.Rectangle((xs[w],ys[w]),20,20,linewidth=2,edgecolor=cols[w],facecolor='none')
        ax1.add_patch(pat)

    for w in range(0,len(cols)):
        pat = patches.Rectangle((xs[w],ys[w]),20,20,linewidth=2,edgecolor=cols[w],facecolor='none')
        ax2.add_patch(pat)
    for w in range(0,len(cols)):
        pat = patches.Rectangle((xs[w],ys[w]),20,20,linewidth=2,edgecolor=cols[w],facecolor='none')
        ax3.add_patch(pat)    
    for w in range(0,len(cols)):
        pat = patches.Rectangle((xs[w],ys[w]),20,20,linewidth=2,edgecolor=cols[w],facecolor='none')
        ax4.add_patch(pat)
    for w in range(0,len(cols)):
        pat = patches.Rectangle((xs[w],ys[w]),20,20,linewidth=2,edgecolor=cols[w],facecolor='none')
        ax5.add_patch(pat)
    for w in range(0,len(cols)):
        pat = patches.Rectangle((xs[w],ys[w]),20,20,linewidth=2,edgecolor=cols[w],facecolor='none')
        ax6.add_patch(pat)

    for i in range(0,len(xs)):
        ax1.text(xs[i]+22,ys[i]+3,stn_in[i], weight = 'bold', fontsize = 20)    


    #tcmap.set_bad('white')
    st = 'Salish Sea Carbonate Chemistry Map, ' + humandate
    plt.suptitle(st,fontsize=20)
    
    fname = rdir + f'{ddmmmyy}_map_' + dss_sig +'.png'

    fig.savefig(fname)
    plt.close()

def surface_buffer_maps(carp,grid, ddmmmyy, rdir,humandate, dss_sig):
    
    #retrieve relevant data for mocsy calculation, calculate mocsy
    tsal = grid.variables['vosaline'][0,0,:,:]
    ttemp = grid.variables['votemper'][0,0,:,:]
    tdic = carp.variables['dissolved_inorganic_carbon'][0,0,:,:]
    tta = carp.variables['total_alkalinity'][0,0,:,:]

    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    tpressure[:] =1
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504
    ttera_is = gsw.t_from_CT(tsra,ttera,tzero)

    response_tup = mocsy.mvars(temp=ttera_is, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tzero, lat=tzero, 
                        optcon='mol/m3', optt='Tinsitu', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    #calculate borate and ohminus concentration
    
    bicarb = hco3
    carb = co3
    #calculate borate, Uppstrom, 1974, looked up in mocsy
    scl=tsra/1.80655
    borat = 0.000232  * scl/10.811
    hplus = 10**(-1*pH)
    borat2 = .0000119*tsra
    ohminus = ttara  - bicarb - 2*carb - borat 

    # - calculates quantities needed for Egleston's factors, and the factors themselves
    
    #Khb is the acidity constant for boric acid - is this an appropriate ref?
    # https://www2.chemistry.msu.edu/courses/cem262/aciddissconst.html
    Khb = 5.81e-10 

    
    S = bicarb + 4*(carb) + (hplus*borat)/(Khb+hplus)+hplus-ohminus
    P = 2*(carb)+bicarb
    AlkC = bicarb + 2*(carb)
    
    DIC = co2+bicarb+carb
    #Alk = bicarb + 2*carb + borat - hplus + ohminus
   
    g_dic = DIC - AlkC**2/S
    b_dic = (DIC*S-AlkC**2)/AlkC
    w_dic = DIC - (AlkC*P)/bicarb
    
    g_alk = (AlkC**2-DIC*S)/AlkC
    b_alk = (AlkC**2/DIC)-S
    w_alk = AlkC - (DIC*bicarb)/P
    
    #### 
    g_dicR = g_dic.reshape(898,398) *1000
    b_dicR = b_dic.reshape(898,398) *1000
    w_dicR = w_dic.reshape(898,398) *-1000
    g_alkR = g_alk.reshape(898,398) *-1000
    b_alkR = b_alk.reshape(898,398) *-1000
    w_alkR = w_alk.reshape(898,398) *1000


    surf_dat = [g_dicR, b_dicR, w_dicR, g_alkR, b_alkR, w_alkR]
    #ranges from nov 13,2014 hindcast. 
    vmins = [-0.7,-0.4,-0.1,-0.4,-0.4,-0.1]
    vmaxs = [0.7,1,0.5,1,1,0.4]
    msk = [1.875e+23,5e+23,6e+23,5e+23,5e+23,2e+23] 

    cl = ['$\gamma_{DIC}$', '$\\beta_{DIC}$', '-$\omega_{DIC}$',\
          '$\gamma_{TA}$', '$\\beta_{TA}$', '-$\omega_{TA}$']
     
    t_cmap = [cm.cm.oxy, cm.cm.oxy, cm.cm.oxy, cm.cm.oxy, cm.cm.oxy, cm.cm.oxy ]
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = \
    plt.subplots(figsize=(20, 27) , nrows=2, ncols=3)

    viz_tools.set_aspect(ax1)
    viz_tools.set_aspect(ax2)
    viz_tools.set_aspect(ax3)
    viz_tools.set_aspect(ax4)
    viz_tools.set_aspect(ax5)
    viz_tools.set_aspect(ax6)

    i = 0
    #'g_dicR',
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,1.875e+23)
    tcmap = t_cmap[i]
    mesh = ax1.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax1)
    cbar.set_label(cl[i], fontsize = 20 )
    ax1.set_title('$CO_{2}$ with DIC',fontsize = 22)
    
    i = 1
    #'b_dicR', 
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,5e+23)
    tcmap = t_cmap[i]
    mesh = ax2.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax2)
    cbar.set_label(cl[i], fontsize = 20 )
    ax2.set_title('pH with DIC',fontsize = 22)

    i = 2
    #'-w_dicR', 
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,6e+23)
    tcmap = t_cmap[i]
    mesh = ax3.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax3)
    cbar.set_label(cl[i], fontsize = 20 )
    ax3.set_title('$\Omega$ with DIC',fontsize = 22)

    i = 3
    #'-g_alkR', 
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,5e+23)
    tcmap = t_cmap[i]
    mesh = ax4.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax4)
    cbar.set_label(cl[i], fontsize = 20 )
    ax4.set_title('$CO_{2}$ with TA',fontsize = 22)
    
    
    i = 4
    #'-b_alkR', 
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,5e+23)
    tcmap = t_cmap[i]
    mesh = ax5.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax5)
    cbar.set_label(cl[i], fontsize = 20 )
    ax5.set_title('pH with TA',fontsize = 22)

    i = 5
    #'w_alkR'
    tplt0 = surf_dat[i]
    tplt = np.ma.masked_values(tplt0,2e+23)
    tcmap = t_cmap[i]
    mesh = ax6.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax6)
    cbar.set_label(cl[i], fontsize = 20 )
    ax6.set_title('$\Omega$ with TA',fontsize = 22)

    #tcmap.set_bad('white')
    st = 'Carbonate Chemistry Buffer Factors, ' + humandate
    plt.suptitle(st,fontsize=20)
    
    fname = rdir + f'{ddmmmyy}_buffmap_' + dss_sig +'.png'

    fig.savefig(fname)
    plt.close()
    
def plume_maps(carp,grid,stns, ddmmmyy, rdir,humandate, dss_sig):
    
    tsal = grid.variables['vosaline'][0,0,:,:]
    ttemp = grid.variables['votemper'][0,0,:,:]
    tdic = carp.variables['dissolved_inorganic_carbon'][0,0,:,:]
    tta = carp.variables['total_alkalinity'][0,0,:,:]

    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    tpressure[:] =1
    tzero = tpressure * 0 
        
    tsra_psu = tsra*35/35.16504

    response_tup = mocsy.mvars(temp=ttera, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tzero, lat=tzero, 
                        optcon='mol/m3', optt='Tpot', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    pHr = pH.reshape(898,398)
    OmA = OmegaA.reshape(898,398)
    
    surf_dat = [tsal, tdic, tta, ttemp, pHr, OmA]
    
    vmins = [15,500,500,5,7.5,0]
    vmaxs = [30,1800,1800,15,8.5,2]
    msk = [0,0,0,0,1e20,1e20]
    cl = ['salinity psu', 'DIC umol/kg', 'TA umol/kg', 'temp deg C', 'pH', 'Omega A']
    t_cmap = [cm.cm.haline, cm.cm.matter, cm.cm.matter, cm.cm.thermal, cm.cm.speed, cm.cm.curl]

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = \
    plt.subplots(figsize=(17, 8.5) , nrows=2, ncols=3)

    viz_tools.set_aspect(ax1)
    viz_tools.set_aspect(ax2)
    viz_tools.set_aspect(ax3)
    viz_tools.set_aspect(ax4)
    viz_tools.set_aspect(ax5)
    viz_tools.set_aspect(ax6)
    
    y1 = 390
    y2 = 460
    x1 = 240
    x2 = 398 
    i = 0
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax1.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax1)
    cbar.set_label(cl[i], fontsize = 20 )
    ax1.set_xticks([])
    ax1.set_yticks([])

    i = 1
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax2.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax2)
    cbar.set_label(cl[i], fontsize = 20 )
    ax2.set_xticks([])
    ax2.set_yticks([])

    i = 2
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax3.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax3)
    cbar.set_label(cl[i], fontsize = 20 )
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    i = 3
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax4.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax4)
    cbar.set_label(cl[i], fontsize = 20 )
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    i = 4
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax5.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax5)
    cbar.set_label(cl[i], fontsize = 20 )
    ax5.set_xticks([])
    ax5.set_yticks([])
    
    i = 5
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax6.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax6)
    cbar.set_label(cl[i], fontsize = 20 )
    ax6.set_xticks([])
    ax6.set_yticks([])
    
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

    #tcmap.set_bad('white')
    st = 'Fraser Plume Carbonate Chemistry, ' + humandate
    plt.suptitle(st,fontsize=20)    
    fname = rdir + f'{ddmmmyy}_plume_' + dss_sig +'.png'

    fig.savefig(fname)
    #plt.show()
    plt.close()    

def parameterspace(pars_profs,stns,rdir,ddmmmyy,humandate, dss_sig):
    fig, ((ax1, ax2, ax3)) = plt.subplots(figsize=(20, 8) , nrows=1, ncols=3)

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
        ax1.legend(stn_in, fontsize = 20)
        ax1.set_title("DIC vs Salinity",fontsize=20)
    #TA
    for i in range(0,nos):
        TA_prof = pars_profs['TA']
        ax2.plot(sal_prof[i,:],TA_prof[i,:],c=cols[i],marker='o',
                           markersize=15, linestyle = '')
        ax2.set_title("TA vs Salinity",fontsize=20)
    #OmegaA
    for i in range(0,nos):
        OmA_prof = pars_profs['OmA']
        ax3.plot(sal_prof[i,:],OmA_prof[i,:],c=cols[i],marker='o',
                           markersize=15, linestyle = '')
        ax3.set_title("OmegaA vs Salinity",fontsize=20)

    ax1.set_xlabel('Salinity (g/kg)', fontsize = 20)
    ax2.set_xlabel('Salinity (g/kg)', fontsize = 20)
    ax3.set_xlabel('Salinity (g/kg)', fontsize = 20)

    ax1.tick_params(axis='y', labelsize=20 )
    ax2.tick_params(axis='y', labelsize=20 )
    ax3.tick_params(axis='y', labelsize=20 )

    ax1.tick_params(axis='x', labelsize=20 )
    ax2.tick_params(axis='x', labelsize=20 )
    ax3.tick_params(axis='x', labelsize=20 )
    
    ax1.set_xlim(15,35)
    ax2.set_xlim(15,35)
    ax3.set_xlim(15,35)
    
    ax1.set_ylim(1800,2300)
    ax2.set_ylim(1800,2300)
    ax3.set_ylim(0,2)
    
    st = 'Salish Sea Carbonate Chemistry Parameterspace Diagram, ' + humandate
    plt.suptitle(st,fontsize=20)
    fname = rdir + f'{ddmmmyy}_pspace_' + dss_sig +'.png'

    fig.savefig(fname)
    plt.close()
    
def days_since_start(date):
    #'13nov14'
    "count days since 01 01 2010, for an easy serial number for plots (01012010 is day 1)"
    daymon = [31,28,31,30,31,30,31,31,30,31,30,31]
    daymon_LY = [31,29,31,30,31,30,31,31,30,31,30,31]
    mon_code = ['jan','feb','mar','apr','may','jun','jul',\
               'aug','sep','oct','nov','dec']

    mon = date[2:5]
    day = int(date[0:2])
    year = int(date[5:7])

    #how many leap days are in the preceding years?
    #if divisible by 4 without remainder, is leapyear
    d, leap = divmod(year, 4)
    leapfactor = d - 3

    #day of year up to & not including first day of month
    index = [i for i, elem in enumerate(mon_code) if mon in elem]
    index = index[0]   
    if leap == 0:
        doy = sum(daymon_LY[0:index])
    if leap != 0:
        doy = sum(daymon[0:index])

    #days in full preceding years    
    yrs_since = year - 10
    dayyear = yrs_since * 365
    
    days_since = dayyear + doy + day + leapfactor

    return days_since

def CCCma_pipe(run_date, prof = True, ncb = True, surfmap = True, \
               buffmap = True, plume = True, pspace = True, pH_T = True):
    #print(run_date)
    #run_date = arrow.get(date)
    ddmmmyy = run_date.format('DDMMMYY').lower()
    humandate = run_date.format('MMM DD, YYYY')
    yyyymmdd = run_date.format('YYYYMMDD')
    
    print('ANALYZING ANALYZING ',humandate)
    print('prof: ', prof, ', ncb: ', ncb, ', surfmap: ', surfmap, 'buffmap: ', buffmap,\
          ', plume: ', plume, ', pspace: ', pspace)
    #change this if you need to change strings
    
    carp1 = f'/results2/SalishSea/hindcast.201812_annex/{ddmmmyy}/SalishSea_1d_{yyyymmdd}_{yyyymmdd}_carp_T.nc'
    grid1 = f'/results2/SalishSea/hindcast.201812_annex/{ddmmmyy}/SalishSea_1d_{yyyymmdd}_{yyyymmdd}_grid_T.nc'
    #carp1 = f'/results/SalishSea/hindcast/{ddmmmyy}/SalishSea_1d_{yyyymmdd}_{yyyymmdd}_carp_T.nc'
    #grid1 = f'/results/SalishSea/hindcast/{ddmmmyy}/SalishSea_1d_{yyyymmdd}_{yyyymmdd}_grid_T.nc'
    carp = nc.Dataset(carp1)
    grid = nc.Dataset(grid1)
    
    rdir = '/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/CCCmaDEV/hind_2018012_feb0119/'
    
    dss = days_since_start(ddmmmyy)
    print(dss)
    dss_sig = str(dss)
    
    pars_profs, stn_list, depths = profiles(carp,grid,cs.STATIONS)
    pars_pts, pt_depths, stn_list2 = point_value(carp,grid,cs.STATIONS)

    if prof == True:
        profile_plotter(pars_profs,depths,cs.STATIONS, humandate, ddmmmyy, rdir, dss_sig)
    if ncb == True:
        ncmaker(stn_list, depths, pt_depths, pars_profs, pars_pts, ddmmmyy, rdir)
    if surfmap == True:
        surface_maps(carp,grid,cs.STATIONS,ddmmmyy,rdir,humandate, dss_sig)
    if buffmap == True:
        surface_buffer_maps(carp,grid,ddmmmyy,rdir,humandate, dss_sig)
    if plume == True:
        plume_maps(carp,grid,cs.STATIONS,ddmmmyy,rdir,humandate, dss_sig)
    if pspace == True:
        parameterspace4(pars_profs,cs.STATIONS,rdir,ddmmmyy,humandate, dss_sig) 
    if plume == True:
        plume_maps(carp,grid,cs.STATIONS,ddmmmyy,rdir,humandate, dss_sig)
    if pH_T == True:
        pH_M = 7
        pH_T = 7.2
        OmA_M = 0.8
        OmA_T = 1
        pCO2_M = 200
        pCO2_T = 400
        pHOmpco2_thres(carp,grid, ddmmmyy, rdir,humandate, dss_sig,pH_M,pCO2_M,OmA_M,pH_T,pCO2_T,OmA_T)
        

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
   

        
def range_analyzer(start,end, prof = True, ncb = True, surfmap = True, \
                   buffmap = True, plume = True, pspace = True, pH_T = True):
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    for i in range(0,len(arrow_array)):
        r = arrow_array[i][0]
        print(r)
        CCCma_pipe(r, prof = prof, ncb = ncb, surfmap = surfmap, \
                   buffmap = buffmap, plume = plume, pspace = pspace, pH_T = pH_T)

        
def pHOmpco2_thres(carp,grid, ddmmmyy, rdir,humandate, dss_sig,pH_M,pCO2_M,OmA_M,pH_T,pCO2_T,OmA_T):
    
    tsal = grid.variables['vosaline'][0,0,:,:]
    ttemp = grid.variables['votemper'][0,0,:,:]
    tdic = carp.variables['dissolved_inorganic_carbon'][0,0,:,:]
    tta = carp.variables['total_alkalinity'][0,0,:,:]

    tsra = np.ravel(tsal)
    ttera = np.ravel(ttemp)
    ttara = np.ravel(tta) * 1e-3
    tdra = np.ravel(tdic) * 1e-3
    tzero = np.zeros_like(tsra)
    tpressure = np.zeros_like(tsra)
    tpressure[:] =1
    tzero = tpressure * 0 

    tsra_psu = tsra*35/35.16504

    response_tup = mocsy.mvars(temp=ttera, sal=tsra_psu, alk=ttara, dic=tdra, 
                       sil=tzero, phos=tzero, patm=tpressure, depth=tzero, lat=tzero, 
                        optcon='mol/m3', optt='Tpot', optp='m',
                        optb = 'l10', optk1k2='m10', optkf = 'dg', optgas = 'Pinsitu')
    pH,pco2,fco2,co2,hco3,co3,OmegaA,OmegaC,BetaD,DENis,p,Tis = response_tup

    pHr = pH.reshape(898,398)
    OmA = OmegaA.reshape(898,398)
    pco2R = pco2.reshape(898,398)

    #surf_dat = [tsal, tdic, tta, ttemp, pHr, OmA]
    surf_dat = [pHr, OmA, pco2R]
    print(OmA[20,250])
    print(pHr[20,250])
    print(pco2R[20,250])

    pH_20 = pH_T - pH_M
    pH_MX = pH_M + (pH_20*5)

    OmA_20 = OmA_T - OmA_M
    OmA_MX = OmA_M + (OmA_20*5)

    pCO2_20 = pCO2_T - pCO2_M
    pCO2_MX = pCO2_M + (pCO2_20*5)

    vmins = [pH_M,OmA_M,pCO2_M]
    vmaxs = [pH_MX,OmA_MX,pCO2_MX]
    msk = [1e20,1e20,1e20]
    cl = ['pH', 'Omega A', 'pCO2']
    t_cmap = [cm.cm.oxy, cm.cm.oxy, cm.cm.oxy]

    fig, ((ax1, ax2, ax3)) = \
    plt.subplots(figsize=(17, 8.5) , nrows=1, ncols=3)

    viz_tools.set_aspect(ax1)
    viz_tools.set_aspect(ax2)
    viz_tools.set_aspect(ax3)


    y1 = 0
    y2 = 898
    x1 = 0
    x2 = 398 
    i = 0
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax1.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax1)
    cbar.set_label(cl[i], fontsize = 20 )
    ax1.set_xticks([])
    ax1.set_yticks([])

    i = 1
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax2.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax2)
    cbar.set_label(cl[i], fontsize = 20 )
    ax2.set_xticks([])
    ax2.set_yticks([])

    i = 2
    tplt0 = surf_dat[i][y1:y2,x1:x2]
    tplt = np.ma.masked_values(tplt0,msk[i])
    tcmap = t_cmap[i]
    mesh = ax3.pcolormesh(tplt, cmap=tcmap, vmin=vmins[i], vmax=vmaxs[i])
    cbar = fig.colorbar(mesh, ax=ax3)
    cbar.set_label(cl[i], fontsize = 20 )
    ax3.set_xticks([])
    ax3.set_yticks([])

    st = 'pH, pCO2, Omega_A threshold plots, ' + humandate
    plt.suptitle(st,fontsize=20)    
    fname = rdir + f'{ddmmmyy}_pHT_' + dss_sig +'.png'

    fig.savefig(fname)
    plt.show()
