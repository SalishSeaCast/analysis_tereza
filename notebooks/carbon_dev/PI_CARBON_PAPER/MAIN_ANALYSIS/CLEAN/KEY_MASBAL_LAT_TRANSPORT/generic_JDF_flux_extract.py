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
    # varname = 'dissolved_inorganic_carbon'

    

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
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*'+ 'dian_U' +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        #print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        if i%50 == 0:
            print(i)
            
    BR_ar = sens_ar
    
    print('done making nclen')
    
    grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
    #[0,:,360:440,20]
    udir = grid['e1t'][0,360:440,20]
    wdir = grid['e3t_0'][0,:,360:440,20]
    print(np.shape(wdir))
    udir_broad = np.zeros([40,80])
    for i in range(0,40):
        udir_broad[i,:] = udir
    tmask = grid['tmask'][0,:,360:440,20]   

    csize_forthis = udir_broad*wdir*tmask  

    for i in range(0,dayslen):

        if (i%50 ==0):
            print(i)

        t_test = nc.Dataset(BR_ar[i])
        #get the UT transport, which is in mmol /s , transfer to mol/day
        tdat = t_test['DIC_UT'][0,:,360:440,20] *(1/1000)*60*60*24
        dic_laytotals = np.nansum(tdat, axis = 1)
        dic_laymeans = np.nansum(tdat, axis = 1)/np.nansum(csize_forthis,axis = 1)        
        
        BR_sums_perday[:,i] =  dic_laytotals
        BR_means_perday[:,i] = dic_laymeans

    fname = sdir_short + '_'+varname+'_JDFfluxsum_perday_alg2.pkl'
    pickle.dump(BR_sums_perday, open(fname, 'wb'))
    fname = sdir_short + '_'+varname+'_JDFfluxmean_perday_alg2.pkl'
    pickle.dump(BR_means_perday, open(fname, 'wb'))
    
    pickle.dump(csize_forthis, open("JDF_transport_cellsize.pkl", 'wb'))
    
    return

    

