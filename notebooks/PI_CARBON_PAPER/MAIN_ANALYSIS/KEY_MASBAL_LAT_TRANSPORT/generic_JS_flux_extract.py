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
    
    
    BR_sums_perday = np.zeros((40,365))
    BR_means_perday = np.zeros((40,365))

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
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*'+ 'dian_V' +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        #print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        if i%50 == 0:
            print(i)
            
    BR_ar = sens_ar
    
    print('done making nclen')
    
    #get 2d dimension (pointing north) of these cells
    grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
    #get udir, vdir, wdir
    vdir = grid['e2t'][0,878,0:120]
    wdir = grid['e3t_0'][0,:,878,0:120]
    print(np.shape(wdir))
    vdir_broad = np.zeros([40,120])
    for i in range(0,40):
        vdir_broad[i,:] = vdir
    tmask = grid['tmask'][0,:,878,0:120]   

    csize_forthis = vdir_broad*wdir*tmask

    for i in range(0,dayslen):

        if (i%50 ==0):
            print(i)

        t_test = nc.Dataset(BR_ar[i])
        #get the UT transport, which is in mmol /s , transfer to mol/day
        tdat = t_test['DIC_VT'][0,:,878,0:120] *(1/1000)*60*60*24
        #print(np.shape(csize_forthis))
        #print(np.shape(tdat))
        dic_laytotals = np.nansum(tdat, axis = 1)
        dic_laymeans = np.nansum(tdat, axis = 1)/np.nansum(csize_forthis,axis = 1)        
        
        BR_sums_perday[:,i] =  dic_laytotals
        BR_means_perday[:,i] = dic_laymeans

    fname = sdir_short + '_'+varname+'_JSfluxsum_perday_alg2.pkl'
    pickle.dump(BR_sums_perday, open(fname, 'wb'))
    fname = sdir_short + '_'+varname+'_JSfluxmean_perday_alg2.pkl'
    pickle.dump(BR_means_perday, open(fname, 'wb'))
    
    pickle.dump(csize_forthis, open("JS_transport_cellsize.pkl", 'wb'))
    
    return

