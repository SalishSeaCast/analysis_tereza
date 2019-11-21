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


ftype = 'ptrc'
ncs_ptrc, arrow_array = make_nclist(start,end,ftype,verbose)

for i in range(0,len(ncs_ptrc)):
    
    tdate = arrow_array[i][0]
    ddmmmyy = tdate.format('DDMMMYY').lower()
    print(ddmmmyy)
    PLACES_withdat = {
    'Deep Bay': {
    'NEMO grid ji': (599, 126),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),


    },

    'Southern Baynes': {
    'NEMO grid ji': (602, 127),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),


    },

    'Northern Baynes': {
    'NEMO grid ji': (646, 127),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Fanny Bay': {
    'NEMO grid ji': (614, 120),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Maple Bay': {
    'NEMO grid ji': (392, 213),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Salt Spring': {
    'NEMO grid ji': (386, 218),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Nanoose Bay': {
    'NEMO grid ji': (517, 190),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Lasqueti Island': {
    'NEMO grid ji': (586, 195),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Main SoG': {
    'NEMO grid ji': (450, 253),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Cortes/Marina': {
    'NEMO grid ji': (732, 157),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Lund/Desolation Sound': {
    'NEMO grid ji': (702, 187),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    'Mouth of Okeover': {
    'NEMO grid ji': (726, 192),
        'ciliates_timeseries':np.zeros([40]),
        'diatoms_timeseries':np.zeros([40]),
        'flagellates_timeseries':np.zeros([40]),

    },

    }

    day = i
    t_ptrc = nc.Dataset(ncs_ptrc[i])



    ciliates = t_ptrc['ciliates']
    diatoms = t_ptrc['diatoms']
    flagellates = t_ptrc['flagellates']

    print('averaging')

    ciliates_daily = np.mean(ciliates,axis=0)
    diatoms_daily = np.mean(diatoms,axis=0)
    flagellates_daily = np.mean(flagellates,axis=0)
    #TA_daily = np.mean(TA,axis=0)
    
    #print('writing to dict')
    for i in range(0,len(list_places)):
        t_place = list_places[i]
        t_j_index = PLACES[t_place]['NEMO grid ji'][0]
        t_i_index = PLACES[t_place]['NEMO grid ji'][1]

        t_ciliates = ciliates_daily[:,t_j_index,t_i_index]
        t_diatoms = diatoms_daily[:,t_j_index,t_i_index]
        t_flagellates = flagellates_daily[:,t_j_index,t_i_index]

        PLACES_withdat[t_place]['ciliates_timeseries'] = t_ciliates
        PLACES_withdat[t_place]['diatoms_timeseries'] = t_diatoms
        PLACES_withdat[t_place]['flagellates_timeseries'] = t_flagellates

        
    str_name = 'EVIE_PLACES_ptrc_'
    t_str = str_name + ddmmmyy + '.nc'
    print('writing to nc')
    print(t_str)

    f = nc.Dataset(t_str, 'w')
    g = f.createGroup('Deep Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Deep Bay']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Deep Bay']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Deep Bay']['flagellates_timeseries'][:]


    g = f.createGroup('Southern Baynes')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Southern Baynes']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Southern Baynes']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Southern Baynes']['flagellates_timeseries'][:]


    g = f.createGroup('Northern Baynes')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Northern Baynes']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Northern Baynes']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Northern Baynes']['flagellates_timeseries'][:]


    g = f.createGroup('Fanny Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Fanny Bay']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Fanny Bay']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Fanny Bay']['flagellates_timeseries'][:]


    g = f.createGroup('Maple Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Maple Bay']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Maple Bay']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Maple Bay']['flagellates_timeseries'][:]

    g = f.createGroup('Salt Spring')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Salt Spring']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Salt Spring']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Salt Spring']['flagellates_timeseries'][:]


    g = f.createGroup('Nanoose Bay')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Nanoose Bay']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Nanoose Bay']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Nanoose Bay']['flagellates_timeseries'][:]

    g = f.createGroup('Lasqueti Island')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lasqueti Island']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lasqueti Island']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lasqueti Island']['flagellates_timeseries'][:]


    g = f.createGroup('Main SoG')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Main SoG']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Main SoG']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Main SoG']['flagellates_timeseries'][:]


    g = f.createGroup('Cortes/Marina')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Cortes/Marina']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Cortes/Marina']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Cortes/Marina']['flagellates_timeseries'][:]


    g = f.createGroup('Lund/Desolation Sound')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lund/Desolation Sound']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lund/Desolation Sound']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Lund/Desolation Sound']['flagellates_timeseries'][:]


    g = f.createGroup('Mouth of Okeover')
    g.createDimension('depths', 40)
    ts = g.createVariable('ciliates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Mouth of Okeover']['ciliates_timeseries'][:]
    ts = g.createVariable('diatoms_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Mouth of Okeover']['diatoms_timeseries'][:]
    ts = g.createVariable('flagellates_timeseries','f4',('depths'))
    ts[:] = PLACES_withdat['Mouth of Okeover']['flagellates_timeseries'][:]



    f.close()

