import arrow
import netCDF4 as nc
import CCCma_profiles as pf
import CCCma_stations as cs
import numpy as np
import gsw
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
from matplotlib import reload


def ncmaker(stn_list, depths, pt_depths, pars_profs, pars_pts, ddmmmyy, rdir):
    
    ncname = rdir + f'{ddmmmyy}_prof.nc'
    f = nc.Dataset(ncname,'w', format='NETCDF4') #'w' stands for write
    g = f.createGroup('Profiles')
    g.createDimension('stn',len(stn_list))
    g.createDimension('z', len(depths))
    g.createDimension('time', None)
    #

    depth = g.createVariable('depth','f4','z')
    depth[:] = depths

    stn = g.createVariable('station','S4','stn')
    for i in range(0,len(stn_list)):
        stn[i] = stn_list[i]
        
    sal_prof = g.createVariable('sal_prof','f4',('stn','z'))
    sal_prof[:] = pars_profs['sal'][:]
    temp_prof = g.createVariable('temp_prof','f4',('stn','z'))
    temp_prof[:] = pars_profs['temp'][:]
    DIC_prof = g.createVariable('DIC_prof','f4',('stn','z'))
    DIC_prof[:] = pars_profs['DIC'][:]
    TA_prof = g.createVariable('TA_prof','f4',('stn','z'))
    TA_prof[:] = pars_profs['TA'][:]
    OmA_prof = g.createVariable('OmA_prof','f4',('stn','z'))
    OmA_prof[:] = pars_profs['OmA'][:]
    pH_prof = g.createVariable('pH_prof','f4',('stn','z'))
    pH_prof[:] = pars_profs['pH'][:]
    O2_prof = g.createVariable('O2_prof','f4',('stn','z'))
    O2_prof[:] = pars_profs['O2'][:]

    sal_profSD = g.createVariable('sal_profSD','f4',('stn','z'))
    sal_profSD[:] = pars_profs['sal_SD'][:]
    temp_profSD = g.createVariable('temp_profSD','f4',('stn','z'))
    temp_profSD[:] = pars_profs['temp_SD'][:]
    DIC_profSD = g.createVariable('DIC_profSD','f4',('stn','z'))
    DIC_profSD[:] = pars_profs['DIC_SD'][:]
    TA_profSD = g.createVariable('TA_profSD','f4',('stn','z'))
    TA_profSD[:] = pars_profs['TA_SD'][:]
    OmA_profSD = g.createVariable('OmA_profSD','f4',('stn','z'))
    OmA_profSD[:] = pars_profs['OmA_SD'][:]
    pH_profSD = g.createVariable('pH_profSD','f4',('stn','z'))
    pH_profSD[:] = pars_profs['pH_SD'][:]
    O2_profSD = g.createVariable('O2_profSD','f4',('stn','z'))
    O2_profSD[:] = pars_profs['O2_SD'][:]

    h = f.createGroup('Station_Averages')
    h.createDimension('stn',len(stn_list))
    h.createDimension('z_pt', len(pt_depths))
    h.createDimension('time', None)

    pt_depth = h.createVariable('depth','f4','z_pt')
    pt_depth[:] = pt_depths

    stn = h.createVariable('station','S4','stn')
    for i in range(0,len(stn_list)):
        stn[i] = stn_list[i]
        
    sal_pt = h.createVariable('sal_pt','f4',('stn','z_pt'))
    sal_pt[:] = pars_pts['sal'][:]
    temp_pt = h.createVariable('temp_pt','f4',('stn','z_pt'))
    temp_pt[:] = pars_pts['temp'][:]
    DIC_pt = h.createVariable('DIC_pt','f4',('stn','z_pt'))
    DIC_pt[:] = pars_pts['DIC'][:]
    TA_pt = h.createVariable('TA_pt','f4',('stn','z_pt'))
    TA_pt[:] = pars_pts['TA'][:]
    OmA_pt = h.createVariable('OmA_pt','f4',('stn','z_pt'))
    OmA_pt[:] = pars_pts['OmA'][:]
    pH_pt = h.createVariable('pH_pt','f4',('stn','z_pt'))
    pH_pt[:] = pars_pts['pH'][:]
    O2_pt = h.createVariable('O2_pt','f4',('stn','z_pt'))
    O2_pt[:] = pars_pts['O2'][:]

    sal_ptSD = h.createVariable('sal_ptSD','f4',('stn','z_pt'))
    sal_ptSD[:] = pars_pts['sal_SD'][:]
    temp_ptSD = h.createVariable('temp_ptSD','f4',('stn','z_pt'))
    temp_ptSD[:] = pars_pts['temp_SD'][:]
    DIC_ptSD = h.createVariable('DIC_ptSD','f4',('stn','z_pt'))
    DIC_ptSD[:] = pars_pts['DIC_SD'][:]
    TA_ptSD = h.createVariable('TA_ptSD','f4',('stn','z_pt'))
    TA_ptSD[:] = pars_pts['TA_SD'][:]
    OmA_ptSD = h.createVariable('OmA_ptSD','f4',('stn','z_pt'))
    OmA_ptSD[:] = pars_pts['OmA_SD'][:]
    pH_ptSD = h.createVariable('pH_ptSD','f4',('stn','z_pt'))
    pH_ptSD[:] = pars_pts['pH_SD'][:]
    O2_ptSD = h.createVariable('O2_ptSD','f4',('stn','z_pt'))
    O2_ptSD[:] = pars_pts['O2_SD'][:]

    f.close()