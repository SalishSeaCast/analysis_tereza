def extractor(start, end, ftype, sdir, varname, fname,  inletmask = False):

    '''USAGE:
    
    for a given variable in raw ncs of the PILA experiment 
    found in directory tree /data/tjarniko/results/BASERUN_EXP
    take by-layer means of the variable throughout the timeperiod specified

    start = '2015-01-01' #start of timeperiod
    end = '2015-12-31' #end of timeperiod (typically a year)
    ftype = 'carp' #type of model result .nc 
    sdir = 'MAIN/BR_3rd_2015' #where under directory tree do we find ncs 
    inletmask = True #are we masking out Toba/Bute/Jervis?
    varname = 'dissolved_inorganic_carbon' #name of variable
    fname = 'BR3_DIC_means' #name of resulting pkl 
    
    import extract_bylayer_mean as ebm
    ebm.extractor(start,end,ftype,sdir, inletmask, varname, fname)
    
    '''
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
    #import gef
    import pickle
        
    #calculate domain size
        
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
    
    #where to store
    daily_means = np.zeros((40,365))
    ncfile_ar = []
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)
    #no of days in array
    dayslen = len(arrow_array)
    
    for i in range(0,dayslen):
        tdate = arrow_array[i][0]
        ddmmmyy = tdate.format('DDMMMYY').lower()
        ymd = tdate.format('YYYYMMDD')
        nc_ncfile = '/data/tjarniko/results/BASERUN_EXP/' + sdir + '/ncs/SKOG_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
        tnc_ncfile = glob.glob(nc_ncfile)
        #print(tnc_ncfile[0])
        ncfile_ar.append(tnc_ncfile[0])
        if i%50 == 0:
            print(i)

    print('done making nclen')

    for i in range(0,dayslen):

        if (i%50 ==0):
            print(ncfile_ar[i])
        t_test = nc.Dataset(ncfile_ar[i])
        
        tdat = t_test[varname][:]
        if i==0:
            print('shape of dataset is')
            print(np.shape(tdat))
        #remove border - 20 grid cells north and west

        tdat[:,:,878:898,:] = np.nan
        tdat[:,:,:,0:20] = np.nan
        #no zeros
        tdat[tdat == 0] = np.nan
        
        if (inletmask == True):
            #mask out inlets
            tdat[0,:,700:898,200:398] = np.nan
            tdat[0,:,550:700,255:398] = np.nan

        if (i==0):
            print('nansum(dataset) -checksum to make sure inlet mask works - false/true should give different answers')
            print(np.nansum(tdat))
        tdat_fc = tdat[0,:,:,:]
        tdat_withvol = tdat_fc*cellsize
        tdat_alldomain = np.nansum(np.nansum(tdat_withvol,axis = 1),axis = 1)


        daily_means[:,i] =  tdat_alldomain

    fname =   fname + '.pkl'
    pickle.dump(daily_means, open(fname, 'wb'))
    
    return