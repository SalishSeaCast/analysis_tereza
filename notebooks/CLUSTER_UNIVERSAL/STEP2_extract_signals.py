import arrow
import netCDF4 as nc
import glob
import numpy as np


## to do this for a different year, change start and end
start ='2012-01-01'
end ='2012-12-31'

#### open the grid file and find cell volumes in m3, and cell thicknesses (for multiplying by mesozooplankton concentrations)
grid = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')

tmask = grid['tmask'][0,:,:,:]
e1t = grid['e1t'][0,:,:]
e2t = grid['e2t'][0,:,:]
e3t_0 = grid['e3t_0'][0,:,:,:]

e1t_proj = np.zeros([40,898,398])
e2t_proj = np.zeros([40,898,398])
for d in range(0,40):
    e1t_proj[d,:] = e1t
    e2t_proj[d,:] = e2t
    
domain_volume = tmask*e1t_proj*e2t_proj*e3t_0

#### Make a list of ncfile-path strings and ymd strings based on a start and end date that you give the function
def get_list_of_model_ncs(start,end, verbose = False, \
                          ncpath = '/results2/SalishSea/nowcast-green.201905',\
                          filetype = 'ptrc_T'):
    "returns a list of model output, need to specify path and filetype (defaults given)"
    
    nclist = []
    ymdlist = []

    start_run = arrow.get(start)
    end_run = arrow.get(end)

    arrow_array = []

    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    dayslen = len(arrow_array)    

    for i in range(0,dayslen):
        if i%50 == 0:
            print(i)
        tdate = arrow_array[i][0]
        ymd = tdate.format('YYYYMMDD')
        ymdlist.append(ymd)
        ncnam = f'{ncpath}/*/SalishSea_1d*{ymd}*{filetype}.nc'            
        t_nc = glob.glob(ncnam)
        nclist.append(t_nc[0])
        
        
    if verbose:
        print(f'first day nc: {nclist[0]}')
        print(f'last day nc: {nclist[-1]}')
        
    return ymdlist, nclist
        
ymdlist, nclist = get_list_of_model_ncs(start,end, verbose = True, \
                          ncpath = '/results2/SalishSea/nowcast-green.201905',\
                         filetype = 'ptrc_T')

### Load the station coordinates found in step 1
stns = nc.Dataset('./DATASETS/X_AND_Y_COORDS.nc')
ycoords = stns['stn_ycoords']
xcoords = stns['stn_xcoords']

### write an extraction function for depth-integrated mesozoo (change for whatever else you want)
def extract_signal_mesozoo(stn_x,stn_y,ncfile):
    
    t_nc = nc.Dataset(ncfile)
    #get depth profile of mesozoo at a given station
    t_mesozoo = (t_nc['mesozooplankton'][0,:,stn_y,stn_x])
    
    #multiply out by cell thickness and sum to get depth_integrated mesozo
    meso_integ = np.nansum(t_mesozoo*e3t_0[:,stn_y,stn_x])
    
    #this returns meso_integ in mmol N / m2
    return meso_integ

### write a looping function called big_extractor that extracts depth-integrated mesozooplankton for each day and station
def big_extractor(nclist, ymdlist,ycoords,xcoords, start_index, end_index, ncstring):
    
    for day in range(start_index,end_index):
        
        #ymd
        t_ymd = ymdlist[day]
        print(t_ymd)
        t_nc = nclist[day]
        
        #for each day, extract the signal for each station
        extracted_signals = np.zeros(len(xcoords))
        
        for stn in range(0,len(xcoords)):
            stn_x = int(xcoords[stn]); stn_y = int(ycoords[stn])
            extracted_signals[stn] = extract_signal_mesozoo(stn_x,stn_y,t_nc)
            
        ### save those signals in an ncfile    
        ncname = f'./DATASETS/{ncstring}_{t_ymd}.nc'
        f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
        #g = f.createGroup('model_output')
        f.createDimension('stn', len(xcoords))
        ts2 = f.createVariable('depthint_mesozoo','f4',('stn'))
        ts2[:] = extracted_signals
        f.close()

## use multiprocessing.Process to run big_extractor in parallel. If 1 day takes 3 minutes and 6 processes are running in parallel, this should take roughly 3 hours?        
ncstring = 'DEPTHINT_MESOZOO'
from multiprocessing import Process
    
def func1():
    print('func1: starting')
    big_extractor(nclist, ymdlist,ycoords,xcoords, 0, 60, ncstring)

def func2():
    print('func2: starting')
    big_extractor(nclist, ymdlist,ycoords,xcoords, 60, 120, ncstring)
      
def func3():
    print('func3: starting')
    big_extractor(nclist, ymdlist,ycoords,xcoords, 120, 180, ncstring)    
    
def func4():
    print('func4: starting')
    big_extractor(nclist, ymdlist,ycoords,xcoords, 180, 240, ncstring)
        
def func5():
    print('func5: starting')
    big_extractor(nclist, ymdlist,ycoords,xcoords, 240, 300, ncstring)
        
def func6():
    print('func6: starting')
    big_extractor(nclist, ymdlist,ycoords,xcoords, 300, 365, ncstring)
    
if __name__ == '__main__':
    p1 = Process(target=func1)
    p1.start()
    p2 = Process(target=func2)
    p2.start()
    p3 = Process(target=func3)
    p3.start()
    p4 = Process(target=func4)
    p4.start()
    p5 = Process(target=func5)
    p5.start()
    p6 = Process(target=func6)
    p6.start()                

    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()
    p6.join()
