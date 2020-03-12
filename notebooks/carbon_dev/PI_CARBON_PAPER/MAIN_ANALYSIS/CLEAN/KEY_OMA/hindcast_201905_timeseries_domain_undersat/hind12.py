yr = '2012'
start = yr + '-01-01'
end = yr + '-12-31'
pkl_name = 'hind_undersat_' + yr + '.pkl'

import pickle
import numpy as np
import netCDF4 as nc
import arrow 

## masking and thresholds

csize = pickle.load(open('../../KEY_MASBAL_LAT_TRANSPORT/pickles/cellsize_alg2.pkl', 'rb'))
mmask = nc.Dataset('/data/tjarniko/MEOPAR/grid/mesh_mask201702.nc')
tmask = mmask['tmask'][0,:,:,:]
#water - wher are cells water?
water = np.where(tmask == 1)
water_dom = (np.shape(tmask[water]))
water_cells = (water_dom[0])

water_vol = np.nansum(np.nansum(np.nansum(csize,axis = 0), axis = 1))
print(water_vol)
print(np.nansum(np.nansum(np.nansum(csize*tmask,axis = 0), axis = 1)))

under_threshold = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]


domain_undersat = np.zeros([11,365])

start_run = arrow.get(start)
end_run = arrow.get(end)

arrow_array = []
for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array.append(r)

dayslen = len(arrow_array)

for i in range(0,dayslen):

    tdate = arrow_array[i][0]
    ymd = tdate.format('YYYYMMDD')

    print(ymd)
    
    t_OmA = nc.Dataset('/data/tjarniko/results/hindcast.201905_dayavg_OmA-pH-pCO2/OmA_plus_'+ymd+'.nc')

    domain_oma = t_OmA['model_output']['OmA'][:,:,:] 

    for t in range(0,len(under_threshold)):
        t_thres = under_threshold[t]
        where_under = np.where((domain_oma<t_thres) & ~(np.isnan(csize)) )
        domain_oma_vol_under = csize[where_under]
        perc_dom_under = np.nansum(domain_oma_vol_under)/water_vol
        domain_undersat[t,i] = perc_dom_under

pickle.dump(domain_undersat, open(pkl_name, 'wb'))
