def generic_extract_fxn(start, end, sdir, sdir_short, varname):

    import matplotlib.pyplot as plt
    import netCDF4 as nc
    import numpy as np
    import scipy as sp
    import datetime as dt
    ""
    from salishsea_tools import (
        nc_tools,
        viz_tools,
        geo_tools,
        tidetools
    )
    import netCDF4 as nc
    import cmocean as cm
    import glob
    import arrow
    import gsw
    import pickle


    sens_ar = []
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    dayslen = len(arrow_array)
    for i in range(0,dayslen):
        tdate = arrow_array[i][0]
        ddmmmyy = tdate.format('DDMMMYY').lower()
        ymd = tdate.format('YYYYMMDD')
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*'+ 'carp_T' +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        #print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        if i%50 == 0:
            print(i)
            
    BR_ar = sens_ar
    
    print('done making nclen')
    
    grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
    #get udir, vdir, wdir
    vdir = grid['e2t'][0,:,:]
    udir = grid['e1t'][0,:,:]
    tmask = grid['tmask'][0,0,:,:]

    surfa = vdir*udir*tmask
    print(np.shape(surfa))
    
    
    
    stor_flx = np.zeros(len(sens_ar))
    mean_flx = np.zeros(len(sens_ar))
    print('calculating air-sea flux')
    i = 0
    for i in range(0,len(sens_ar)):
        if i%50 == 0:
            print(i)
        G = nc.Dataset(sens_ar[i])
        var_tmp = G.variables['co2_flux_mmol_m2_s'][0,0:878,20:398]
        var_tmp[var_tmp == 1e+20] = 0
        var_tmp2 = var_tmp * surfa[0:878,20:398]
        totflux = np.sum(np.sum(var_tmp2))
        meanflux = np.nanmean(var_tmp)
        #mol/day
        totflux_daily_moles = totflux * 60 * 60 * 24 * (1/1000)
        meanflux_daily_moles = meanflux * 60 * 60 * 24 * (1/1000)
        stor_flx[i] = totflux_daily_moles
        mean_flx[i] = meanflux_daily_moles




    fname = sdir_short + '_ASflux_sum_perday_alg2.pkl'
    pickle.dump(stor_flx, open(fname, 'wb'))
    fname = sdir_short + '_ASflux_mean_perday_alg2.pkl'
    pickle.dump(mean_flx, open(fname, 'wb'))
    
    pickle.dump(surfa[0:878,20:398], open("ASflux_cellsize.pkl", 'wb'))
    
    return

