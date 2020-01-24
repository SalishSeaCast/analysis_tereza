def generic_extract_fxn(start, end, ftype, sdir, sdir_short, varname):

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
    import gef
    import pickle
    
    


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
        nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
        tnc_sens = glob.glob(nc_sens)
        #print(tnc_sens[0])
        sens_ar.append(tnc_sens[0])
        if i%50 == 0:
            print(i)
            
    BR_ar = sens_ar
    
    print('done making nclen')

    for i in range(0,365):

        if (i%5 ==0):
            print(i)

        t_test = nc.Dataset(BR_ar[i])
        tdat = t_test[varname][:]
        tdat[:,878:898,:] = np.nan
        tdat[:,:,0:20] = np.nan
        tdat[tdat == 0] = np.nan
        tdat_fc = tdat[0,:,:,:]
        tdat_alldomain = np.zeros([40])
        for q in range(0,40):
            tdat_alldomain[q] = np.nanmean(tdat_fc[q,:,:])
        
        
        if (i%5 ==0):
            print('fixed tdat_alldomain???')
            print('assigned nans, they really should show up')
            print(tdat_alldomain)
        BR_means_perday[:,i] =  tdat_alldomain


    fname = sdir_short + '_'+varname+'_means_perday_alg2.pkl'
    pickle.dump(BR_means_perday, open(fname, 'wb'))

    
    return

    