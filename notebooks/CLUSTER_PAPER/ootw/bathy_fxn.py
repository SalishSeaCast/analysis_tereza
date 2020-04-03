'''POINT: 
| CALLS: 
| NOTES: 
| TAKES : 
| GIVES: 
| USAGE: 
| USED IN: 
| written 201*, tjšj '''

def threshold_value(bathy, threshold):
    '''given a bathymetry file, returns a numpy array of size SSNEMO of 0s (land) 1s (shallower than threshold, and 2s (deeper than threshold)'''
    import numpy as np 
    
    depth = np.array(bathy.variables['Bathymetry'])
    depth_bin = depth
    walrus = np.where(np.logical_and(depth>0, depth<threshold))
    seal = np.where(depth >= threshold)
    depth_bin[walrus] = 1
    depth_bin[seal] =2

    return depth_bin

def stn_threshold_list(bathy, d_stn_y, d_stn_x, threshold):
    '''gives a list telling, for each station, whether it's deeper (2) or shallower(1) than threshold'''
    import numpy as np
    
    status = np.zeros(len(d_stn_x))
    for i in range(0,len(d_stn_x)):
        
        depth_bin = threshold_value(bathy, threshold)
        if depth_bin[d_stn_y[i],d_stn_x[i]] ==2 :
            status[i] = 2
        if depth_bin[d_stn_y[i],d_stn_x[i]] ==1 :
            status[i] = 1     
        
    return status

def closest_deep_station(this_y,this_x,depths):
    '''SUPERCEDED by closest_deep_station_withtrace
    finds the closest deep station,
    given the x and y indices of a shallow station
    '''
    import numpy as np
    
    ys, xs  = np.where(depths ==2)
    dists = []
    for i in range(0,len(xs)):
        dist = np.sqrt(abs(this_x-xs[i])**2+(this_y-ys[i])**2)
        dists.append(dist)
    closest = min(dists)
    ind =dists.index(min(dists))
    #ind = ind[0]
    x_close = xs[ind]
    y_close = ys[ind]

    return x_close, y_close, closest

def master_stn_list(bathy,stn_x,stn_y,threshold):
    '''given a bathymetry .nc file, list of x and y indices of stations, and a threshold depth, 
    returns:
    stn_status, which is 1 if the stn is shallower than the given depth, 2 if deep enough
    close_x, which is 999 if the station is deep enough, or lists x index of the closest deep station
    close_y, same logic
    and distance (in x-y index space) to nearest stn, which is unimportant but useful for checking'''
    import numpy as np
    
    close_x = np.zeros_like(stn_x, dtype= float)
    close_y = np.zeros_like(stn_x, dtype = float)
    dist_to_nearest_deep_stn = np.zeros_like(stn_x, dtype = float)
    close_x[:] = 999
    close_y[:] = 999
    dist_to_nearest_deep_stn[:] = 999
    
    stn_status = stn_threshold_list(bathy, stn_y, stn_x, threshold)
    where_shallow = np.where(stn_status == 1)
    depth_bin = threshold_value(bathy,threshold)
    where_shallow = np.squeeze(where_shallow)
    for i in range(0,len(where_shallow)):
        #print(i)
        ind = where_shallow[i]
        x, y, closest = closest_deep_station(stn_y[ind],stn_x[ind],depth_bin)
        close_x[ind] = x
        close_y[ind] = y
        dist_to_nearest_deep_stn[ind] = closest
        
    return stn_status, close_x, close_y, dist_to_nearest_deep_stn

def find_closest_cell_to_given_depth2(grid,depth):
    '''for a given bathymetry file and given depth, 
    finds the z-index of the cell closest to that depth
    superceded, uses different grid?'''
    
    import numpy as np
    
    delt_z = grid.e3t
    z_spac = (delt_z[:,:,0,0])
    z_spac = np.squeeze(z_spac)
    depth_spac = np.zeros_like(z_spac, dtype = float)
    for i in range(0,len(depth_spac)):
        depth_spac[i] = sum(z_spac[0:i])
    
    finder = abs(depth_spac - depth)
    close = np.min(finder)
    
    ind = np.where(finder == close)
    ind = np.squeeze(ind)
    return ind


def find_closest_cell_to_given_depth(trc,depth):
    '''for a given tracer dataset path and given depth, finds the z-index of the cell closest to that depth '''
    import numpy as np
    import netCDF4 as nc
    trd = nc.Dataset(trc)
    depthar = np.array(trd.variables['deptht'])
    #print(depthar)
    depth_spac = np.zeros_like(depth, dtype = float)
    finder = abs(depth-depthar)
    close = np.min(finder)
    
    ind = np.where(finder == close)
    ind = np.squeeze(ind)
    return ind

def closest_deep_station_withtrace(this_y,this_x,depths, ncDATA):
    '''finds the closest deep station, given the x and y indices of a shallow station
    checks with model output array to make sure '''
    import numpy as np
    import netCDF4 as nc
    
    #retreive salinity
    sal = np.array(ncDATA.variables['vosaline'])
    sal = np.squeeze(sal)
    surf_sal = sal[0,:,:]
    
    #are the datasets the same size?
    walrus = (surf_sal.shape == depths.shape)
    print(walrus)
    
    #eliminate where surface salinity is zero/no data
    ys, xs  = np.where((depths ==2) & (surf_sal != 0))
    dists = []
    for i in range(0,len(xs)):
        dist = np.sqrt(abs(this_x-xs[i])**2+(this_y-ys[i])**2)
        dists.append(dist)
    closest = min(dists)
    ind =dists.index(min(dists))

    x_close = xs[ind]
    y_close = ys[ind]

    
    border = 0 
    if (x_close ==0) :
        border =1
        print('yr at the border & you have a problem')
    if (y_close ==0) :
        border =1
        
    return x_close, y_close, closest

def closest_deep_station_with_array(this_y,this_x,depths, sal):
    '''finds the closest deep station, given the x and y indices of a shallow station
    checks with model output array to make sure '''
    import numpy as np

    surf_sal = sal[0,:,:]
    
    #are the datasets the same size?
    walrus = (surf_sal.shape == depths.shape)
    #print(walrus)
    
    #eliminate where surface salinity is zero/no data
    ys, xs  = np.where((depths ==2) & (surf_sal != 0))
    dists = []
    for i in range(0,len(xs)):
        dist = np.sqrt(abs(this_x-xs[i])**2+(this_y-ys[i])**2)
        dists.append(dist)
    closest = min(dists)
    ind =dists.index(min(dists))

    x_close = xs[ind]
    y_close = ys[ind]

    
    border = 0 
    if (x_close ==0) :
        border =1
        print('yr at the border & you have a problem')
    if (y_close ==0) :
        border =1
        
    return x_close, y_close, closest