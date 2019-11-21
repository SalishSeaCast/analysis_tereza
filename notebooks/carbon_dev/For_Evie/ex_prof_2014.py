import netCDF4 as nc
import numpy as np
import arrow
import glob

### change this for a different year
start = '2014-01-01'
end = '2014-12-31'
verbose = False
### 

list_places = ['Deep Bay', 'Southern Baynes', 'Northern Baynes', \
               'Fanny Bay', 'Maple Bay', 'Salt Spring', 'Nanoose Bay',\
               'Lasqueti Island', 'Main SoG', 'Cortes/Marina', \
               'Lund/Desolation Sound', 'Mouth of Okeover']

PLACES = {
'Deep Bay': {

# deg E, deg N

'lon lat': (-124.7392, 49.4606),

# indices of nearest NEMO model grid point

# j is the latitude (y) direction, i is the longitude (x) direction

'NEMO grid ji': (599, 126),

},

'Southern Baynes': {

'lon lat': (-124.7457, 49.4760),

'NEMO grid ji': (602, 127),

},

'Northern Baynes': {

'lon lat': (-124.8924, 49.6492),

'NEMO grid ji': (646, 127),

},

'Fanny Bay': {

'lon lat': (-124.8227, 49.5086),

'NEMO grid ji': (614, 120),

},

'Maple Bay': {

'lon lat': (-123.5947, 48.8140),

'NEMO grid ji': (392, 213),

},

'Salt Spring': {

'lon lat': (-123.5513, 48.7993),

'NEMO grid ji': (386, 218),

},

'Nanoose Bay': {

'lon lat': (-124.1359, 49.2609),

'NEMO grid ji': (517, 190),

},

'Lasqueti Island': {

# deg E, deg N

'lon lat': (-124.3384, 49.5442),

'NEMO grid ji': (586, 195),

},

'Main SoG': {

'lon lat': (-123.5832, 49.1177),

'NEMO grid ji': (450, 253),

},

'Cortes/Marina': {

'lon lat': (-125.0194, 50.0418),

'NEMO grid ji': (732, 157),

},

'Lund/Desolation Sound': {

'lon lat': (-124.7666, 49.9804),

'NEMO grid ji': (702, 187),

},

'Mouth of Okeover': {

'lon lat': (-124.8174, 50.0805),

'NEMO grid ji': (726, 192),

},

}


def make_nclist(start,end,ftype,verbose):

    print('Variable type is '+ftype)
    print('')
    nc_array = []
    start_run = arrow.get(start)
    end_run = arrow.get(end)
    
    
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    dayslen = len(arrow_array)
    
    for i in range(0,dayslen):
        tdate = arrow_array[i][0]
        if verbose == True:
            print('DATE IS: ' +str(tdate))
        ddmmmyy = tdate.format('DDMMMYY').lower()
        if verbose == True:
            print('Arrow ddmmmyy format is: '+ str(ddmmmyy))
        ymd = tdate.format('YYYYMMDD')
        if verbose == True:
            print('Arrow ymd format is: '+ str(ymd))
            
        t_nc = '/results/SalishSea/hindcast.201905/' \
        + str(ddmmmyy)+ '/SalishSea_*' + str(ymd) + '_' + str(ymd) + '_' + ftype + '_T.nc'
        t_ncname = glob.glob(t_nc)
        if verbose == True:
            print('file is: ' + t_ncname[0])
            print('')
        nc_array.append(t_ncname[0])
        
    return nc_array, arrow_array


ftype = 'carp'
ncs_carp, arrow_array = make_nclist(start,end,ftype,verbose)

ftype = 'grid'
ncs_grid, arrow_array = make_nclist(start,end,ftype,verbose)

for i in range(0,len(ncs_grid)):
    
    tdate = arrow_array[i][0]
    ddmmmyy = tdate.format('DDMMMYY').lower()
    print(ddmmmyy)
    PLACES_withdat = {
    'Deep Bay': {
    'NEMO grid ji': (599, 126),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),

    },

    'Southern Baynes': {
    'NEMO grid ji': (602, 127),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),

    },

    'Northern Baynes': {
    'NEMO grid ji': (646, 127),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Fanny Bay': {
    'NEMO grid ji': (614, 120),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Maple Bay': {
    'NEMO grid ji': (392, 213),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Salt Spring': {
    'NEMO grid ji': (386, 218),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Nanoose Bay': {
    'NEMO grid ji': (517, 190),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Lasqueti Island': {
    'NEMO grid ji': (586, 195),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Main SoG': {
    'NEMO grid ji': (450, 253),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Cortes/Marina': {
    'NEMO grid ji': (732, 157),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Lund/Desolation Sound': {
    'NEMO grid ji': (702, 187),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    'Mouth of Okeover': {
    'NEMO grid ji': (726, 192),
        'T_timeseries':np.zeros([40]),
        'S_timeseries':np.zeros([40]),
        'DIC_timeseries':np.zeros([40]),
        'TA_timeseries':np.zeros([40]),
    },

    }

    day = i
    t_carp = nc.Dataset(ncs_carp[i])
    t_grid = nc.Dataset(ncs_grid[i])

    SAL = t_grid['vosaline']
    TEMP = t_grid['votemper']
    DIC = t_carp['dissolved_inorganic_carbon']
    TA = t_carp['total_alkalinity']

    print('averaging')

    SAL_daily = np.mean(SAL,axis=0)
    TEMP_daily = np.mean(TEMP,axis=0)
    DIC_daily = np.mean(DIC,axis=0)
    TA_daily = np.mean(TA,axis=0)
    print('writing to dict')
    for i in range(0,len(list_places)):
        t_place = list_places[i]
        t_j_index = PLACES[t_place]['NEMO grid ji'][0]
        t_i_index = PLACES[t_place]['NEMO grid ji'][1]

        t_DIC = DIC_daily[:,t_j_index,t_i_index]
        t_TA = TA_daily[:,t_j_index,t_i_index]
        t_SAL = SAL_daily[:,t_j_index,t_i_index]
        t_TEMP = TEMP_daily[:,t_j_index,t_i_index]

        PLACES_withdat[t_place]['T_timeseries'] = t_TEMP
        PLACES_withdat[t_place]['S_timeseries'] = t_SAL
        PLACES_withdat[t_place]['DIC_timeseries'] = t_DIC
        PLACES_withdat[t_place]['TA_timeseries'] = t_TA
        
    str_name = 'EVIE_PLACES_'
    t_str = str_name + ddmmmyy + '.nc'
    print('writing to nc')
    print(t_str)

    f = nc.Dataset(t_str, 'w')
    g = f.createGroup('Deep Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Deep Bay']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Deep Bay']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Deep Bay']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Deep Bay']['TA_timeseries'][:]

    g = f.createGroup('Southern Baynes')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Southern Baynes']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Southern Baynes']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Southern Baynes']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Southern Baynes']['TA_timeseries'][:]

    g = f.createGroup('Northern Baynes')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Northern Baynes']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Northern Baynes']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Northern Baynes']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Northern Baynes']['TA_timeseries'][:]

    g = f.createGroup('Fanny Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Fanny Bay']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Fanny Bay']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Fanny Bay']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Fanny Bay']['TA_timeseries'][:]

    g = f.createGroup('Maple Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Maple Bay']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Maple Bay']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Maple Bay']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Maple Bay']['TA_timeseries'][:]

    g = f.createGroup('Salt Spring')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Salt Spring']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Salt Spring']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Salt Spring']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Salt Spring']['TA_timeseries'][:]

    g = f.createGroup('Nanoose Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Nanoose Bay']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Nanoose Bay']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Nanoose Bay']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Nanoose Bay']['TA_timeseries'][:]

    g = f.createGroup('Lasqueti Island')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lasqueti Island']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lasqueti Island']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lasqueti Island']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lasqueti Island']['TA_timeseries'][:]

    g = f.createGroup('Main SoG')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Main SoG']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Main SoG']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Main SoG']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Main SoG']['TA_timeseries'][:]

    g = f.createGroup('Cortes/Marina')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Cortes/Marina']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Cortes/Marina']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Cortes/Marina']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Cortes/Marina']['TA_timeseries'][:]

    g = f.createGroup('Lund/Desolation Sound')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lund/Desolation Sound']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lund/Desolation Sound']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lund/Desolation Sound']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lund/Desolation Sound']['TA_timeseries'][:]

    g = f.createGroup('Mouth of Okeover')
    g.createDimension('depths', 40)
    ts = g.createVariable('T_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Mouth of Okeover']['T_timeseries'][:]
    ts = g.createVariable('S_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Mouth of Okeover']['S_timeseries'][:]
    ts = g.createVariable('DIC_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Mouth of Okeover']['DIC_timeseries'][:]
    ts = g.createVariable('TA_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Mouth of Okeover']['TA_timeseries'][:]


    f.close()
