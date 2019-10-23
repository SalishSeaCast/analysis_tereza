import arrow
import netCDF4 as nc
import glob

import CCCma_profiles as pf
import CCCma_stations as cs
import CCCma_ncmaker as nm
import CCCma_pspace as ps
import numpy as np
import gsw
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import reload
import cmocean as cm

from salishsea_tools import (
#    nc_tools,
    viz_tools
#    geo_tools,
#    tidetools
)

def range_analyzer(tdir, start,end, rdir, prof = True, ncb = True, surfmap = True, \
                   buffmap = True, plume = True, pspace = True, pH_T = True):
    
    """Wrapper for CCCma plotting code
    TJSJ, MOAD group, UBC, February 2019
    
    Input keywords:
      start,end, prof = True, ncb = True, surfmap = True, \
      buffmap = True, plume = True, pspace = True, pH_T = True):
    
    Uses arrow to give start, end dates and call CCCma_pipe, 
    which calls individual plotting codes depending on user input
        
    Option to make 6 plots: 
        1) prof: depth profiles of rel. conc. 
        2) surfmap: surface concentration map
        3) buffmap: map of buffer factors (doesn't work right now)
        4) plume: surface concentration in the plume region
        5) pspace: set of 4 parameterspace plots
        6) pH_T: threshold plots for pH, OmA, pCO2

    Option to make a netcdf file of relevant quantities and their std deviations:
        1) ncb: (short for netcdf builder) 
    
    Usage: 
    
    start = '2017-05-01'
    end = '2017-09-01'

    CCCma.range_analyzer(start,end, prof = True, ncb = True, surfmap = True, \
                   buffmap = True, plume = True, pspace = True, pH_T = True)

    """

    start_run = arrow.get(start)
    end_run = arrow.get(end)
    arrow_array = []
    for r in arrow.Arrow.span_range('day', start_run, end_run):
        arrow_array.append(r)

    for i in range(0,len(arrow_array)):
        r = arrow_array[i][0]
        print(r)
        CCCma_pipe(tdir, r, rdir, prof = prof, ncb = ncb, surfmap = surfmap, \
                   buffmap = buffmap, plume = plume, pspace = pspace, pH_T = pH_T)
        
def days_since_start(date):
    #'13nov14'
    "count days since 01 01 2010, for an easy serial number for plots (01012010 is day 1)"
    daymon = [31,28,31,30,31,30,31,31,30,31,30,31]
    daymon_LY = [31,29,31,30,31,30,31,31,30,31,30,31]
    mon_code = ['jan','feb','mar','apr','may','jun','jul',\
               'aug','sep','oct','nov','dec']

    mon = date[2:5]
    day = int(date[0:2])
    year = int(date[5:7])

    #how many leap days are in the preceding years?
    #if divisible by 4 without remainder, is leapyear
    d, leap = divmod(year, 4)
    leapfactor = d - 3

    #day of year up to & not including first day of month
    index = [i for i, elem in enumerate(mon_code) if mon in elem]
    index = index[0]   
    if leap == 0:
        doy = sum(daymon_LY[0:index])
    if leap != 0:
        doy = sum(daymon[0:index])

    #days in full preceding years    
    yrs_since = year - 10
    dayyear = yrs_since * 365
    
    days_since = dayyear + doy + day + leapfactor

    return days_since
        
        
def CCCma_pipe(tdir, run_date, rdir, prof = True, ncb = True, surfmap = True, \
               buffmap = True, plume = True, pspace = True, pH_T = True):
    #print(run_date)
    #run_date = arrow.get(date)
    ddmmmyy = run_date.format('DDMMMYY').lower()
    humandate = run_date.format('MMM DD, YYYY')
    yyyymmdd = run_date.format('YYYYMMDD')
    
    print('ANALYZING ANALYZING ',humandate)
    print('prof: ', prof, ', ncb: ', ncb, ', surfmap: ', surfmap, 'buffmap: ', buffmap,\
          ', plume: ', plume, ', pspace: ', pspace)
    #change this if you need to change strings
    
    carp1 = f'/{tdir}/SKOG_1d*_carp_T_{yyyymmdd}-{yyyymmdd}.nc'
    carp2 = glob.glob(carp1)
    grid1 = f'/{tdir}/SKOG_1d*_grid_T_{yyyymmdd}-{yyyymmdd}.nc'
    grid2 = glob.glob(grid1)
    #print(grid1)
    #print(grid2)
    #print(carp1)
    #print(carp2)
    carp = nc.Dataset(carp2[0])
    grid = nc.Dataset(grid2[0])
    
    
    dss = days_since_start(ddmmmyy)
    print(dss)
    dss_sig = str(dss)
    
    print('walrus if statements')
    if ((prof == True) | (ncb == True) | (pspace == True)):
        reload(pf)

        pars_pts, pt_depths, stn_list2 = pf.point_value(carp,grid,cs.STATIONS)
        pars_profs, stn_list, depths = pf.profiles(carp,grid,cs.STATIONS)
        
        if ncb == True:
            nm.ncmaker(stn_list, depths, pt_depths, pars_profs, pars_pts, ddmmmyy, rdir)

        if pspace == True:
            reload(ps)

            ps.parameterspace4(pars_profs,cs.STATIONS,rdir,ddmmmyy,humandate, dss_sig) 

        if prof == True:
            print('prof')
            pf.profile_plotter(pars_profs,depths,cs.STATIONS, humandate, ddmmmyy, rdir, dss_sig)
    #
    if surfmap == True:
        reload(mp)
        mp.surface_maps(carp,grid,cs.STATIONS,ddmmmyy,rdir,humandate, dss_sig)
    #if buffmap == True:
    #    surface_buffer_maps(carp,grid,ddmmmyy,rdir,humandate, dss_sig)
    if plume == True:
        mp.plume_maps(carp,grid,cs.STATIONS,ddmmmyy,rdir,humandate, dss_sig)

    #if pH_T == True:
    #    pH_M = 7
    #    pH_T = 7.2
    #    OmA_M = 0.8
    #    OmA_T = 1
    #    pCO2_M = 200
    #    pCO2_T = 400
    #    pHOmpco2_thres(carp,grid, ddmmmyy, rdir,humandate, dss_sig,pH_M,pCO2_M,OmA_M,pH_T,pCO2_T,OmA_T)