# Copyright 2013-2016 The Salish Sea MEOPAR contributors
# and The University of British Columbia

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#: Places to extract for Evie 
PLACES = {
    # Baynes Sound
    'Deep Bay': {
        # deg E, deg N
        'lon lat': (-124.7392, 49.4606),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (599, 126),
    },
    'Southern Baynes': {
        # deg E, deg N
        'lon lat': (-124.7457, 49.4760),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (602, 127),
    },
    'Northern Baynes': {
        # deg E, deg N
        'lon lat': (-124.8924, 49.6492),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (646, 127),
    },
    'Fanny Bay': {
        # deg E, deg N
        'lon lat': (-124.8227, 49.5086),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (614, 120),
    },
    # South SoG
    'Maple Bay': {
        # deg E, deg N
        'lon lat': (-123.5947, 48.8140),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (392, 213),
    },
    'Salt Spring': {
        # deg E, deg N
        'lon lat': (-123.5513, 48.7993),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (386, 218),
    },
    'Nanoose Bay': {
        # deg E, deg N
        'lon lat': (-124.1359, 49.2609),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (517, 190),
    },
    'Lasqueti Island': {
        # deg E, deg N
        'lon lat': (-124.3384, 49.5442),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (586, 195),
    },
    'Main SoG': {
        # deg E, deg N
        'lon lat': (-123.5832, 49.1177),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (450, 253),
    },
    # Discovery Islands
    'Cortes/Marina': {
        # deg E, deg N
        'lon lat': (-125.0194, 50.0418),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (732, 157),
    },
    'Lund/Desolation Sound': {
        # deg E, deg N
        'lon lat': (-124.7666, 49.9804),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (702, 187),
    },
    'Mouth of Okeover': {
        # deg E, deg N
        'lon lat': (-124.8174, 50.0805),
        # indices of nearest NEMO model grid point
        # j is the latitude (y) direction, i is the longitude (x) direction
        'NEMO grid ji': (726, 192),
    },

}

def DispGeoLocs():
    """Display locations in map coordinates

    :returns: figure handle
    :rtype: :py:class:`matplotlib.figure.Figure`
    """
    from mpl_toolkits.basemap import Basemap
    from matplotlib import pyplot as plt
    places2=PLACES.copy()
    width = 160000; lon_0 = -124.5; lat_0 = 49.4
    fig=plt.figure(figsize=(12,12))
    m = Basemap(width=width,height=width,projection='aeqd', resolution='h',
                lat_0=lat_0,lon_0=lon_0)
    m.drawmapboundary()
    m.drawcoastlines(linewidth=0.5)
    m.drawrivers()
    m.drawparallels(range(40,60,1))
    m.drawmeridians(range(-130,-110,1))
    #plt.title('EC River Stations')

    # map stations:
    for pl in places2.keys():
        if 'lon lat' in places2[pl].keys():
            lon,lat=places2[pl]['lon lat']
            if (47<lat<51) & (-128<lon<-120):
                if pl in ('Sandy Cove','Calamity Point','Port Moody','Vancouver','New Westminster','Delta DDL node','East node','Boundary Bay','Duke Pt.'):
                    xpt, ypt = m(lon, lat)
                    xpt2, ypt2 = m(lon+.03, lat)
                    m.plot(xpt,ypt,'ro')
                    plt.text(xpt2,ypt2,pl,fontsize=10,fontweight='bold',
                            ha='left',va='center',color='r')
                else:
                    xpt, ypt = m(lon, lat)
                    xpt2, ypt2 = m(lon-.03, lat)
                    m.plot(xpt,ypt,'ro')
                    plt.text(xpt2,ypt2,pl,fontsize=10,fontweight='bold',
                            ha='right',va='center',color='r')
    return fig,m


def DispGridLocs(mesh_mask='/ocean/eolson/MEOPAR/NEMO-forcing/grid/mesh_mask201702_noLPE.nc'):
    """Display locations in NEMO model grid coordinates

    :arg mesh_mask: string with path to the meshmask you would like to plot; 201702 default
    :type mesh_mask: str

    :returns: figure handle
    :rtype: :py:class:`matplotlib.figure.Figure`
    """
    import numpy as np # only grid
    from mpl_toolkits.basemap import Basemap
    import netCDF4 as nc # only grid
    from matplotlib import pyplot as plt
    from salishsea_tools import viz_tools # only grid
    places2=PLACES.copy()
    with nc.Dataset(mesh_mask) as fm:
        tmask=np.copy(fm.variables['tmask'])
        e3t_0=np.copy(fm.variables['e3t_0'])
    bathy=np.sum(e3t_0[0,:,:,:]*tmask[0,:,:,:],0)
    cm=plt.cm.get_cmap('Blues')
    cm.set_bad('lightgray')
    fig,ax=plt.subplots(1,2,figsize=(18,18))
    viz_tools.set_aspect(ax[0])
    viz_tools.set_aspect(ax[1])
    ax[0].pcolormesh(np.ma.masked_where(tmask[0,0,:,:]==0,bathy),cmap=plt.cm.get_cmap('Blues'))
    # map stations:
    for pl in places2.keys():
        if 'NEMO grid ji' in places2[pl].keys() and places2[pl]['NEMO grid ji'] is not None:
            j,i=places2[pl]['NEMO grid ji']
            if pl in ('Sandy Cove','Calamity Point','Port Moody','Vancouver','New Westminster',
                      'East Node','Duke Pt.','Halibut Bank','Cherry Point','Central SJDF','Friday Harbor'):
                ax[0].plot(i,j,'ro')
                ax[0].text(i+4,j,pl,fontsize=10,fontweight='bold',
                        ha='left',va='center',color='r')
            elif pl in ('Sand Heads','Delta DDL node','Central Node','Delta BBL node','Cluster_9',
                      'East Node','Woodwards Landing'):
                ax[0].plot(i,j,'ro')
            else:
                ax[0].plot(i,j,'ro')
                ax[0].text(i-4,j,pl,fontsize=10,fontweight='bold',
                        ha='right',va='center',color='r')
    xl=(240,340)
    yl=(400,450)
    ax[1].pcolormesh(np.ma.masked_where(tmask[0,0,:,:]==0,bathy),cmap=plt.cm.get_cmap('Blues'))
    # map stations:
    for pl in places2.keys():
        if 'NEMO grid ji' in places2[pl].keys() and places2[pl]['NEMO grid ji'] is not None:
            j,i=places2[pl]['NEMO grid ji']
            if (xl[0]<i<xl[1]) & (yl[0]<j<yl[1]):
                if pl in ('Sandy Cove','Calamity Point','Port Moody','Vancouver','New Westminster','Friday Harbor',
                          'East Node','Duke Point','Halibut Bank','Cherry Point', 'Sand Heads','Cluster_9','Central SJDF'):
                    ax[1].plot(i,j,'ro')
                    ax[1].text(i+1,j,pl,fontsize=10,fontweight='bold',
                            ha='left',va='center',color='r')
                elif pl in ('Delta BBL node'):
                    ax[1].plot(i,j,'ro')
                    ax[1].text(i,j-2,pl,fontsize=10,fontweight='bold',
                            ha='center',va='center',color='r')
                elif pl in ('Delta DDL node'):
                    ax[1].plot(i,j,'ro')
                    ax[1].text(i,j+2,pl,fontsize=10,fontweight='bold',
                            ha='right',va='center',color='r')
                else:
                    ax[1].plot(i,j,'ro')
                    ax[1].text(i-1,j,pl,fontsize=10,fontweight='bold',
                            ha='right',va='center',color='r')
    ax[1].set_xlim(xl[0],xl[1])
    ax[1].set_ylim(yl[0],yl[1])
    return fig


