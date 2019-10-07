import xarray as xr
import glob
import numpy as np

def import_bathy(filepath):
    filepath = '/data/tjarniko/MEOPAR/NEMO-forcing/grid/mesh_mask_SalishSea2.nc'
    grid = xr.open_dataset(filepath)
    return grid

def make_stns(spacing):
    xs =[]
    ys =[]
    count = 0
    stn_x = []
    stn_y = []
    for i in range (0,200):
        x = 2+i*spacing
        if(x>397):
            break
        else:
            xs.append(x)
        
    for j in range(0,400):
        y = 2+j*spacing
        if(y>897):
            break
        else:
            ys.append(y)
            
    for a in range(0,len(xs)):
        for b in range(0,len(ys)):
            ts_x = xs[a]
            ts_y = ys[b]
            stn_x.append(ts_x)
            stn_y.append(ts_y)
            
    return stn_x, stn_y

def filter_stn_in_domain(stn_x,stn_y,fmask):
    d_stn_x = []
    d_stn_y = []
    for s in range(0,len(stn_x)):
        x = stn_x[s]
        y = stn_y[s]
        stn = fmask.values[y-1:y+2,x-1:x+2]
        #if there are no zeroes in s
        if((0 in stn) == False):
            d_stn_x.append(x)
            d_stn_y.append(y)
            
    return d_stn_x, d_stn_y

def create_tracelist(monlist,daylist,nday,year):
    '''monlist = ['jan','feb','mar','apr','may','jun','jul','aug','sep',
      'oct','nov','dec']
      daylist = [31,28,31,30,31,30,31,31,30,31,30,31]'''
    if (year == 2013):
        ye = 13
    if (year == 2014):
        ye = 14
    if (year == 2016):
        ye = 16
    if (year == 2015):
        ye = 15
    if (year == 2017):
        ye = 17
    if (year == 2018):
        ye = 18
    #should write an exception here)
    doy = []
    for i in range(0,12):
        mon = monlist[i]
        day = daylist[i]
        for j in range(1,day+1):
            dom = j
            if (dom<10):
                dom = str(0) + str(dom)
            else:
                dom = str(dom)
            daypath = dom + monlist[i] + str(ye)
            doy.append(daypath)

    trace_list = []
    
        #make paths to dataset results for first 310 days of the year
    for i in range(0,nday):
        Day_OI = doy[i]
        
        #print('The day of interest is ' + Day_OI)
        #
        #respath = '/results/SalishSea/hindcast'
        #print(year)
        respath = 'walrus'
        #print(i)
        respath = '/results/SalishSea/spinup.201905/'
        print(respath)
#        if ((year == 2015) | ((year == 2016) & (i < 325))) :
#            respath = '/results/SalishSea/hindcast.201812'
#        if ((year == 2016) & (i > 324)) :
#            respath = '/results2/SalishSea/hindcast.201812_annex'
#        if ((year == 2017) | (year == 2018)) :    
#            respath = '/results2/SalishSea/hindcast.201812_annex'
        #print(respath)
        daypath = respath + '/' + Day_OI 
        #print('xx')
        print(daypath)
        tracers = 'SalishSea_1d*ptrc_T.nc'
        ptt = daypath +'/' + tracers
        trace = glob.glob(ptt)
        #print(trace)
        #print(i)
        #print('***')
        #print(trace)

        trace = (trace[0])
        trace_list.append(trace)
    
    return trace_list

def create_physlist(monlist,daylist,nday,year):
    
    if (year == 2013):
        ye = 13
    if (year == 2014):
        ye = 14
    if (year == 2016):
        ye = 16
    if (year == 2015):
        ye = 15
    if (year == 2017):
        ye = 17
    if (year == 2018):
        ye = 18
    doy = []
    respath = '/results/SalishSea/hindcast'
    for i in range(0,12):
        mon = monlist[i]
        day = daylist[i]
        for j in range(1,day+1):
            dom = j
            if (dom<10):
                dom = str(0) + str(dom)
            else:
                dom = str(dom)
            daypath = dom + monlist[i] + str(ye)
            doy.append(daypath)

    trace_list = []
    
        #make paths to dataset results for first 310 days of the year
    for i in range(0,nday):
        Day_OI = doy[i]
        #print('The day of interest is ' + Day_OI)
        #respath = '/results/SalishSea/hindcast'
        respath = 'walrus'
        respath = '/results/SalishSea/spinup.201905/'
        #print(i)
#        if ((year == 2015) | ((year == 2016) & (i < 325))) :
#           respath = '/results/SalishSea/hindcast.201812'
#        if ((year == 2016) & (i > 324)) :
#            respath = '/results2/SalishSea/hindcast.201812_annex'
#        if ((year == 2017) | (year == 2018)) :    
#            respath = '/results2/SalishSea/hindcast.201812_annex'
        daypath = respath + '/' + Day_OI 
        #print(daypath)
        tracers = 'SalishSea_1d*grid_T.nc'
        ptt = daypath +'/' + tracers
        trace = glob.glob(ptt)
        trace = (trace[0])
        #print(trace)
        trace_list.append(trace)
    
    return trace_list

def create_physlist_W(monlist,daylist,nday,year):
    doy = []
    if (year == 2013):
        ye = 13
    if (year == 2014):
        ye = 14
    if (year == 2016):
        ye = 16
    if (year == 2015):
        ye = 15
    if (year == 2017):
        ye = 17
    if (year == 2018):
        ye = 18
    
    print(ye)
    doy = []
    respath = '/results/SalishSea/hindcast'
    for i in range(0,12):
        mon = monlist[i]
        day = daylist[i]
        for j in range(1,day+1):
            dom = j
            if (dom<10):
                dom = str(0) + str(dom)
            else:
                dom = str(dom)
            daypath = dom + monlist[i] + str(ye)
            doy.append(daypath)

    trace_list = []
    
    for i in range(0,nday):
        Day_OI = doy[i]
        #print('The day of interest is ' + Day_OI)
        #respath = '/results/SalishSea/hindcast'
        respath = 'walrus'
        #print(i)
        respath = '/results/SalishSea/spinup.201905/'
        daypath = respath + '/' + Day_OI 
        #print(daypath)
        tracers = 'SalishSea_1d*grid_W.nc'
        ptt = daypath +'/' + tracers
        trace = glob.glob(ptt)
        trace = (trace[0])
        #print(trace)
        trace_list.append(trace)
    
    return trace_list


