def yearnit_de(spacing,stn_b,stn_e,year):
    
    import map_fxn as mf
    import time
    import xarray as xr
    import numpy as np
    print('Spacing between stations: ' + str(spacing))
    
    bath = '/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])
    
    stn_x, stn_y = mf.make_stns(spacing)
    d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)
    
    no_stns = len(d_stn_x)
    
    monlist = ['jan','feb','mar','apr','may','jun','jul','aug','sep',
          'oct','nov','dec']
    if year == 2016:
        daylist = [31,29,31,30,31,30,31,31,30,31,30,31]
        noday = 366
    else:
        daylist = [31,28,31,30,31,30,31,31,30,31,30,31]
        noday = 365

    # a list of all the model outputs (tracers) for this year
    trace_list = mf.create_tracelist(monlist,daylist,noday,year)
    
    print("Number of stations:" + str(no_stns))
    
    for stn in range(stn_b,stn_e):
    
        print('station is: ' ,str(stn))
        print('x is :', d_stn_x[stn])
        print('y is :', d_stn_y[stn])
        
        ts_x = d_stn_x[stn]
        ts_y = d_stn_y[stn]
        
        daily_nit = np.zeros(noday)
        
        end_day = noday
        for day in range(1,end_day+1):
            if (day%20 ==0) :
                print(day)
            trc = trace_list[day-1]
            #print(trc)
            nemo = xr.open_mfdataset(trc)
            #print(nemo)
            no3 = np.array(nemo['nitrate'])
            no3 = np.squeeze(no3)
            surf_no3 = no3[0:3,ts_y,ts_x]
            surf_no3 = sum(surf_no3)
            daily_nit[day-1] = surf_no3
            
        nit = xr.Dataset({'surface_nitrogen':(['t'], daily_nit)})
        dirname = '/ocean/tjarniko/MEOPAR/NC_hindcast/' + str(year) +'/NIT_TS/stn_'
        stn_name = dirname + str(stn)  + 'surf_nit_sp' + str(spacing)+ '.nc'
        nit.to_netcdf(stn_name)

        
        