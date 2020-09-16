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
        nc_ncfile = '/data/tjarniko/results/RIV_PIL/' + sdir + '/SalishSeaCast_1d_*'+ ftype +'*' + ymd + '-' + ymd + '.nc'
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
        tdat_alldomain = np.zeros([40])
        for q in range(0,40):
            tdat_alldomain[q] = np.nanmean(tdat_fc[q,:,:])

        daily_means[:,i] =  tdat_alldomain

    fname =  fname + '.pkl'
    pickle.dump(daily_means, open(fname, 'wb'))
    
    return