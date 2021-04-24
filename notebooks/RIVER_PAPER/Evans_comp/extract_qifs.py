import arrow
import pickle
import numpy as np
import netCDF4 as nc
import glob

qifs_j = 761;
qifs_i = 137
### loop for extraction of results
start ='2015-01-01'
end ='2020-12-31'

hrs = np.arange(0,24,1)/24
#print(hrs)
hr_seg_of_yr = hrs/366

start_run = arrow.get(start)
end_run = arrow.get(end)

arrow_array = []

for r in arrow.Arrow.span_range('day', start_run, end_run):
    arrow_array.append(r)
print('done')
dayslen = len(arrow_array)

date_doy = np.zeros([dayslen*24])
temp_doy = np.zeros([dayslen*24])
sal_doy = np.zeros([dayslen*24])
DIC_doy = np.zeros([dayslen*24])
TA_doy = np.zeros([dayslen*24])

count = 0

for i in range(0,dayslen):

    tdate = arrow_array[i][0]
    ymd = tdate.format('YYYYMMDD')
    doy = tdate.format('DDD')
    yy = tdate.format('YYYY')
    mm = tdate.format('MMM').lower()
    sy = tdate.format('YY')
    dd = tdate.format('DD')
    dayfolder = f'{dd}{mm}{sy}'

    if (( int(yy) == 2016) | ( int(yy) == 2020)):
        daysinyear = 366
    else:
        daysinyear = 365
    print(dayfolder)
    
    w = f'/results2/SalishSea/nowcast-green.201905/{dayfolder}/SalishSea_1h_{ymd}_{ymd}_grid_T.nc'
    w2 = glob.glob(w)
    gridT = nc.Dataset(w2[0])
    temp = gridT['votemper'][:,0,qifs_j,qifs_i]
    sal = gridT['vosaline'][:,0,qifs_j,qifs_i]
    temp_doy[count:count+24] = temp
    sal_doy[count:count+24] = sal
    
    w = f'/results2/SalishSea/nowcast-green.201905/{dayfolder}/SalishSea_1h_{ymd}_{ymd}_carp_T.nc'
    w2 = glob.glob(w)
    carpT = nc.Dataset(w2[0])
    DIC = carpT['dissolved_inorganic_carbon'][:,0,qifs_j,qifs_i]
    TA = carpT['total_alkalinity'][:,0,qifs_j,qifs_i]
    DIC_doy[count:count+24] = DIC
    TA_doy[count:count+24] = TA 
    
    
    datehrs = int(yy)+(int(doy)-1)/daysinyear + hr_seg_of_yr
    date_doy[count:count+24] = datehrs

    count = count+24
    
pickle.dump(date_doy, open("./pkls/date_doy.pkl", 'wb'))
pickle.dump(DIC_doy, open("./pkls/DIC_doy.pkl", 'wb'))
pickle.dump(TA_doy, open("./pkls/TA_doy.pkl", 'wb'))
pickle.dump(temp_doy, open("./pkls/temp_doy.pkl", 'wb'))
pickle.dump(sal_doy, open("./pkls/sal_doy.pkl", 'wb'))