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
    
    
    grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
    #get udir, vdir, wdir
    vdir = grid['e2t'][0,:,:]
    udir = grid['e1t'][0,:,:]
    wdir = grid['e3t_0'][0,:,:,:]
    tmask = grid['tmask'][0,:,:,:]

    # w = np.array([[2,3],[2,3]])
    # x = np.array([[2,4],[2,3]])
    # print(w*x)
    surfa = vdir*udir
    surfa_broad = np.zeros([40,898,398])
    for i in range(0,40):
        surfa_broad[i,:,:] = surfa

    csize_recalc = surfa_broad*wdir*tmask
    csize_recalc[:,878:898,:] = 0
    csize_recalc[:,:,0:20] = 0
    csize_recalc[csize_recalc==0] = np.nan
    cellsize = csize_recalc


    BR_sums_perday = np.zeros((40,365))
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
        tdat_withvol = tdat_fc*cellsize
        tdat_alldomain = np.nansum(np.nansum(tdat_withvol,axis = 1),axis = 1)
        if (i%5 ==0):
            print(tdat_withvol[:,250,250])
            print('assigned nans, they really should show up')
            print(tdat_alldomain)
        BR_sums_perday[:,i] =  tdat_alldomain


    fname = sdir_short + '_'+varname+'_sums_perday_alg2.pkl'
    pickle.dump(BR_sums_perday, open(fname, 'wb'))
    pickle.dump(cellsize, open("cellsize_alg2.pkl", 'wb'))
    
    return

    
# def make_nclen(start,end,ftype, sdir):
#     base_ar = []
#     sens_ar = []
#     start_run = arrow.get(start)
#     end_run = arrow.get(end)
#     arrow_array = []
#     for r in arrow.Arrow.span_range('day', start_run, end_run):
#         arrow_array.append(r)

#     dayslen = len(arrow_array)
#     for i in range(0,dayslen):
#         tdate = arrow_array[i][0]
#         ddmmmyy = tdate.format('DDMMMYY').lower()
#         ymd = tdate.format('YYYYMMDD')
#         nc_sens = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
#         tnc_sens = glob.glob(nc_sens)
#         #print(tnc_sens[0])
#         sens_ar.append(tnc_sens[0])
        
#     return sens_ar
# start = '2015-01-01'
# end = '2015-12-31'
# ftype = 'carp'
# sdir = 'MAIN/BR_2nd_2015'
# sdir_short = 'BR_2nd_2015'
# varname = 'dissolved_inorganic_carbon'
