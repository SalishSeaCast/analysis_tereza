import arrow
import netCDF4 as nc
import CCCma_profiles as pf
import CCCma_stations as cs
import CCCma_ncmaker as nm
import CCCma_pspace as ps
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
    
    #b = 0
    for s in stns:
        stn_list.append(s)
        
        
        #print('Calculating point values for ' + stns[s]['fullname'])
        ty = stns[s]['y']
        tx = stns[s]['x']
        b = stns[s]['serialno']

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

            tsr[tsr == 0] = np.nan
            tter[tter == 0] = np.nan
            tdr[tdr == 0] = np.nan
            ttar[ttar == 0] = np.nan

            sal_pt[b,d] = np.nanmean(tsr, axis = 1)
            sal_ptSD[b,d] = np.nanstd(tsr, axis = 1)

            temp_pt[b,d] = np.nanmean(tter, axis = 1)
            temp_ptSD[b,d] = np.nanstd(tter, axis = 1)

            DIC_pt[b,d] = np.nanmean(tdr, axis = 1)
            DIC_ptSD[b,d] = np.nanstd(tdr, axis = 1)

            TA_pt[b,d] = np.nanmean(ttar, axis = 1)
            TA_ptSD[b,d] = np.nanstd(ttar, axis = 1)

            O2_pt[b,d] = np.nanmean(to2r, axis = 1)
            O2_ptSD[b,d] = np.nanstd(to2r, axis = 1)

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
        
            pH_pt[b,d] = np.nanmean(pHra, axis = 1)
            pH_ptSD[b,d] = np.nanstd(pHra, axis = 1)

            OmA_pt[b,d] = np.nanmean(OmegaAra, axis = 1)
            OmA_ptSD[b,d] = np.nanstd(OmegaAra, axis = 1)
            
            sal_pt[sal_pt == 0 ] = np.nan
            temp_pt[temp_pt == 0 ] = np.nan
            DIC_pt[DIC_pt == 0 ] = np.nan
            TA_pt[TA_pt == 0 ] = np.nan
            pH_pt[pH_pt == 0 ] = np.nan
            OmA_pt[OmA_pt == 0 ] = np.nan

    
    pars_pts = {'sal': sal_pt, 'sal_SD': sal_ptSD,\
    'temp': temp_pt, 'temp_SD': temp_ptSD,\
    'DIC': DIC_pt, 'DIC_SD' : DIC_ptSD,\
    'TA': TA_pt, 'TA_SD' : TA_ptSD,\
    'OmA': OmA_pt, 'OmA_SD': OmA_ptSD,\
    'pH': pH_pt, 'pH_SD': pH_ptSD,\
    'O2' : O2_pt, 'O2_SD': O2_ptSD}
    
    return pars_pts, pt_depths, stn_list 

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
    
    for s in stns:
        stn_list.append(s)
        stn = stns[s]
        #print('Calculating profiles for ' + stns[s]['fullname'])
        ty = stns[s]['y']
        tx = stns[s]['x']
        b = stns[s]['serialno']

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

        tsr[tsr == 0] = np.nan
        tter[tter == 0] = np.nan
        tdr[tdr == 0] = np.nan
        ttar[ttar == 0] = np.nan

        sal_prof[b,:] = np.nanmean(tsr, axis = 1)
        sal_profSD[b,:] = np.nanstd(tsr, axis = 1)
        
        temp_prof[b,:] = np.nanmean(tter, axis = 1)
        temp_profSD[b,:] = np.nanstd(tter, axis = 1)

        DIC_prof[b,:] = np.nanmean(tdr, axis = 1)
        DIC_profSD[b,:] = np.nanstd(tdr, axis = 1)

        TA_prof[b,:] = np.nanmean(ttar, axis = 1)
        TA_profSD[b,:] = np.nanstd(ttar, axis = 1)

        O2_prof[b,:] = np.nanmean(to2r, axis = 1)
        O2_profSD[b,:] = np.nanstd(to2r, axis = 1)

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

        pH_prof[b,:] = np.nanmean(pHra, axis = 1)
        pH_profSD[b,:] = np.nanstd(pHra, axis = 1)

        OmA_prof[b,:] = np.nanmean(OmegaAra, axis = 1)
        OmA_profSD[b,:] = np.nanstd(OmegaAra, axis = 1)
        
        #depth = w.variables['deptht'][:]
        sal_prof[sal_prof == 0 ] = np.nan
        temp_prof[temp_prof == 0 ] = np.nan
        DIC_prof[DIC_prof == 0 ] = np.nan
        TA_prof[TA_prof == 0 ] = np.nan
        pH_prof[pH_prof == 0 ] = np.nan
        OmA_prof[OmA_prof == 0 ] = np.nan
        
    
        
    pars_profs = {'sal': sal_prof, 'sal_SD': sal_profSD,\
         'temp': temp_prof, 'temp_SD': temp_profSD,\
         'DIC': DIC_prof, 'DIC_SD' : DIC_profSD,\
         'TA': TA_prof, 'TA_SD' : TA_profSD,\
         'OmA': OmA_prof, 'OmA_SD': OmA_profSD,\
         'pH': pH_prof, 'pH_SD': pH_profSD,\
         'O2' : O2_prof, 'O2_SD': O2_profSD}
    return pars_profs, stn_list, prof_depth

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
