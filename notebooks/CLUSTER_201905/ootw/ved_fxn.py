def yearved_de(spacing,stn_b,stn_e,year):
    
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
        print(noday)
    if year == 2015:
        daylist = [31,28,31,30,31,30,31,31,30,31,30,31]
        noday = 365
        print(noday)
        

    # a list of all the model outputs (tracers) for this year    
    trace_list = mf.create_physlist_W(monlist,daylist,noday,year)
    
    print("Number of stations:" + str(no_stns))
    
    for stn in range(stn_b,stn_e):
    
        print('station is: ' ,str(stn))
        print('x is :', d_stn_x[stn])
        print('y is :', d_stn_y[stn])
        
        ts_x = d_stn_x[stn]
        ts_y = d_stn_y[stn]
        
        daily_ved = np.zeros(noday)
        
        end_day = noday
        for day in range(1,end_day+1):
            if (day%20 ==0) :
                print(day)
            trc = trace_list[day-1]
            #print(trc)
            nemo = xr.open_mfdataset(trc)
            #print(nemo)
            ved = np.array(nemo.variables['vert_eddy_diff'])
            ved = np.squeeze(ved)
            
            depth = np.array(nemo.variables['depthw_bounds'])
            depth = depth[:,1]
            
            ved_here = ved[:,ts_y,ts_x]
            ved_here = np.squeeze(ved_here)
            where0 = np.where(ved_here == 0)
            where0 = np.squeeze(where0)
            first_oob = (where0[1])#surface is also always zero
            depth_here = depth[first_oob-1]
            total_ved = np.sum(ved_here[0:first_oob])
            avg_ved = total_ved/depth_here
            
            daily_ved[day-1] = avg_ved
    
        ved = xr.Dataset({'daily_ved':(['t'], daily_ved)})
        stn_name = '/ocean/tjarniko/MEOPAR/NC_hindcast/' + str(year) + '/VED_TS/stn_' + str(stn)  + 'avg_ved_sp' + str(spacing)+ '.nc'
        ved.to_netcdf(stn_name)
