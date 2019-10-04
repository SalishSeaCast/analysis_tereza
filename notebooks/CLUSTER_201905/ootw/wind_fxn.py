'''POINT: | CALLS: | NOTES: | TAKES : | GIVES: | USAGE: | USED IN: | written 201*, tjšj '''

def create_GEMlist(year):   
    ''' POINT: creates a python list of atmospheric forcing files for a given year | CALLS: | NOTES: path to forcing files is hardcoded (currently GEM2.5 operational), can be changed with dirstring. Doesn't believe in leap years. Will only work with years for which all atmospheric forcing files are available | TAKES: a year int | GIVES: a list of strings specifying forcing .nc files for every day in that year, eg 'ops_y2016m01d01.nc' | USAGE: gemlist = create_GEMlist(2016) | ROLE IN CLUSTER PROJECT: | written 2017, tjšj ''' 
    daylist = [31,28,31,30,31,30,31,31,30,31,30,31]
    if year == 2016 :
        daylist = [31,29,31,30,31,30,31,31,30,31,30,31]
    gemlist = []
   
    print(year) 
    for i in range(0,12):
        dirstring = '/results/forcing/atmospheric/GEM2.5/operational/'
        filestart = dirstring + 'ops_y' + str(year) + 'm'
        mon = i+1
        
        if mon < 10:
            monstr = str(0) + str(mon)
        else:
            monstr = str(mon)
        
        filestart = filestart+monstr + 'd'
        #print(filestart)
        
        day = daylist[i]
        for j in range(1,day+1):
            dom = j
            if (dom<10):
                dom = str(0) + str(dom)
            else:
                dom = str(dom)
            #print(dom)
            gem = filestart + dom + '.nc'
            gemlist.append(gem)
    
    print('this should be the number of days')
    print(len(gemlist))
    return gemlist


def avg_daily_windmag_windstress_energy(windfile):
    '''POINT: calculates average daily wind magnitudes, wind stresses (wind magnitudes squared) and wind energies (wind magnitudes cubed), from a given atmospheric forcing file |CALLS: | NOTES: | TAKES : a string specifying the location of an .nc GEM2.5 operational atm forcing file for a single day, which has dimension: 24 hrs * 266 lats * 256 lons | GIVES: 3 arrays of avg daily winds, windstresses, wind energies dimensions 266 lats * 256 lons | USAGE: daily_wm, daily_ws, daily_we  = avg_daily_windmags(windfile) | ROLE IN CLUSTER PROJECT: | written 2017, tjšj '''
    import xarray as xr
    import numpy as np
    wind = xr.open_dataset(windfile)

    wind_lats = wind.nav_lat
    wind_lons = wind.nav_lon
    u_wind = wind.u_wind
    v_wind = wind.v_wind
    
    wind_mag = np.sqrt((u_wind)*(u_wind)+(v_wind)*(v_wind))
    wind_stress = wind_mag ** 2
    wind_energy = wind_mag ** 3
    
    daily_wm = np.average(wind_mag, axis=0)
    daily_ws = np.average(wind_stress, axis=0)
    daily_we = np.average(wind_energy, axis=0)
    return daily_wm, daily_ws, daily_we

def get_atmogrid():
    '''POINT: returns lats and lons of an atmospheric forcing grid | CALLS: | NOTES: uses a sample .nc file as an example  | TAKES : nothing | GIVES: two meshgrids, wind_lats, wind_lons | USAGE: wind_lats, wind_lons = get_atmogrid() | ROLE IN CLUSTER PROJECT: | written 2017, tjšj '''
    import xarray as xr
    wind = xr.open_dataset('/results/forcing/atmospheric/GEM2.5/operational/ops_y2016m01d01.nc')
    wind_lats = wind.nav_lat
    wind_lons = wind.nav_lon
    return wind_lats, wind_lons

def get_oceangrid():
    '''POINT: returns lats and lons of the model salish sea configuration grid | CALLS: | NOTES: uses the grid meshmask .nc file as an example  | TAKES : nothing | GIVES: two meshgrids, ocean_lats, ocean_lons | USAGE: ocean_lats, ocean_lons = get_oceangrid() | ROLE IN CLUSTER PROJECT: | written 2017, tjšj '''
    import xarray as xr
    ocean = xr.open_dataset('/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc')
    ocean_lats = ocean.nav_lat
    ocean_lons = ocean.nav_lon
    return ocean_lats, ocean_lons

def weights_intrp_mtx(odata):
    '''POINT: interpolates data from the atmospheric grid onto the ocean grid using a weights matrix | CALLS: | NOTES: this is much faster than interpolating with scipy! | TAKES : one array of size lats 266 x lons 256 of atmospheric data | GIVES: an array of NEMO-Salish size of the new data interpolated | USAGE: | ROLE IN CLUSTER PROJECT: | written by Michael Dunphy, adapted by tjšj, 2017'''
    from IPython import embed
    import netCDF4 as nc
    import scipy.interpolate as spi
    import scipy.sparse as sp
    import numpy as np
    '''takes a matrix of size lon(256), lat(266) and interpolates it onnto the NEMO gird
    size is 266 x 256'''
    weightsfile = '/home/mdunphy/MEOPAR/NEMO-forcing/grid/weights-gem2.5-ops.nc'
    with nc.Dataset(weightsfile) as f:
        s1 = f.variables['src01'][:]-1  # minus one for fortran-to-python indexing
        s2 = f.variables['src02'][:]-1
        s3 = f.variables['src03'][:]-1
        s4 = f.variables['src04'][:]-1
        w1 = f.variables['wgt01'][:]
        w2 = f.variables['wgt02'][:]
        w3 = f.variables['wgt03'][:]
        w4 = f.variables['wgt04'][:]
       
    NO = odata.size   # number of operational grid points
    NN = s1.size      # number of NEMO grid points
    
    # Build matrix
    n = np.array([x for x in range(0,NN)])
    M1 = sp.csr_matrix((w1.flatten(), (n, s1.flatten())), (NN,NO))
    M2 = sp.csr_matrix((w2.flatten(), (n, s2.flatten())), (NN,NO))
    M3 = sp.csr_matrix((w3.flatten(), (n, s3.flatten())), (NN,NO))
    M4 = sp.csr_matrix((w4.flatten(), (n, s4.flatten())), (NN,NO))
    M = M1+M2+M3+M4
    
        # Interpolate by matrix multiply - quite fast
    ndata = M*odata.flatten()

    # Reshape to NEMO shaped array
    ndata=ndata.reshape(s1.shape)
    
    return ndata

def produce_interpolated_wind_netcdf(year, daystart, dayend):
    '''POINT: produces a netcdf file containing SalishSeaNEMO sized arrays of wind magnitude, wind stress, and wind energy for each in a given series of days | CALLS: create_gemlist, weights_intrp_mtx (both from wind_fxn) | NOTES: day numbering starts at 0! the fname is hardcoded to a given directory into which nc files are deposited, could be changed maybe | TAKES : year, day start (starts at 0), day end  | GIVES: netcdfs in a given directory| USAGE: produce_interpolated_wind_netcdf | ROLE IN CLUSTER PROJECT: | written 2017, tjšj '''
    
    import xarray as xr
    import map_fxn as mf
    
    gemlist = create_GEMlist(year)
    for i in range(daystart,dayend+1):
        if (i%10 == 0):
            print(i)
        f = gemlist[i]
        daily_wm, daily_ws, daily_we = avg_daily_windmag_windstress_energy(f)
        daily_wm_int = weights_intrp_mtx(daily_wm)
        daily_ws_int = weights_intrp_mtx(daily_ws)
        daily_we_int = weights_intrp_mtx(daily_we)
        
        
        windthings = xr.Dataset({'daily_avg_windmag': (['y','x'], daily_wm_int),
                                 'daily_avg_windstress': (['y','x'], daily_ws_int),
                                 'daily_avg_windenergy': (['y','x'], daily_we_int),})
        
        fname = '/ocean/tjarniko/MEOPAR/WINDFILES_interp/YEAR' + str(year) + '_day' +str(i) +'.nc' 
        
        windthings.to_netcdf(fname)
        
        
def yearwinds_de(spacing,stn_b,stn_e,year):
    '''POINT: for a given set of stations, produces 3 timeseries of wind forcing data (magnitude, stress, energy)  station as lists in a single netcdf file (one netcdf file per station) | CALLS: import_bathy, make_stns, filster_station_in_domain (from map_fxn), reads in netcdf files produced by produce_interpolated_wind_netcdf | NOTES: days start at 0, stations start at 0, for a spacing of 10 there should be 580 stations, the start and end station option is so that it can be run in parallel if necessary, directory into which nc files are deposited is hardcoded under stn_name | TAKES : spacing (has been decided as 10 for this project), year, start station, end station | GIVES: netcdf files that contain 3 windrelated timeseries for each station, in a given repo | USAGE: yearwinds_de(10,1,580,2016) | ROLE IN CLUSTER PROJECT: | written 2017, tjšj '''
    
    import map_fxn as mf
    import time
    import xarray as xr
    import numpy as np
    print('Spacing between stations: ' + str(spacing))
    
    noday = 365
    if year == 2016:
        noday = 366
    
    print('no. days, from yearwinds_de fxn')
    print(noday)
    
    bath = '/results/nowcast-sys/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = mf.import_bathy(bath)
    fmask = (grid.fmask[0,0,:,:])
    
    stn_x, stn_y = mf.make_stns(spacing)
    d_stn_x, d_stn_y = mf.filter_stn_in_domain(stn_x,stn_y,fmask)
    print('number of stns: ' + str(len(d_stn_x)))
    
    for stn in range(stn_b,stn_e):
    
        print('station is: ' ,str(stn))
        print('x is :', d_stn_x[stn])
        print('y is :', d_stn_y[stn])
        
        ts_x = d_stn_x[stn]
        ts_y = d_stn_y[stn]
        
        daily_windmag = np.zeros(noday)
        daily_windstress = np.zeros(noday)
        daily_windenergy = np.zeros(noday)
        
        end_day = noday
        
        
        
        for day in range(0,end_day):
            ##open the winds file
            #dayfile
            
            dirstring = '/ocean/tjarniko/MEOPAR/WINDFILES_interp/'
            #fname = '/ocean/tjarniko/MEOPAR/WINDFILES_interp/'
            filestart = dirstring + 'YEAR' + str(year) + '_day' + str(day) + '.nc'
           
            #print(day)
            #print(filestart) 
            wf = xr.open_dataset(filestart)
            ## extract wm, ws, we for the station for the day
            ## attach it to the corresponding array above
            ## save the resulting 3 lists to a netcdf file as below
            
            #avg_daily_windmag_windstress_energy
            
            daily_we = wf.daily_avg_windenergy[ts_y,ts_x]
            daily_we = daily_we.values
            daily_ws = wf.daily_avg_windstress[ts_y,ts_x]
            daily_ws = daily_ws.values
            daily_wm = wf.daily_avg_windmag[ts_y,ts_x]
            daily_wm = daily_wm.values
            ##be gentle about OB1 errors!!!!
            

            
            daily_windmag[day] = daily_wm
            daily_windstress[day] = daily_ws
            daily_windenergy[day] = daily_we
            
            
            winds = xr.Dataset({'wind_mags':(['t'], daily_windmag), 'wind_stresses':(['t'], daily_windstress), 
                                'wind_energy':(['t'], daily_windenergy)})
            stn_name = '/ocean/tjarniko/MEOPAR/NC_hindcast/' + str(year) + '/WIND_TS/stn_' + str(stn) + '_wind_data_sp' + str(spacing) + '.nc' 
            winds.to_netcdf(stn_name)
            
            
            
            
            
            
