import netCDF4 as nc
import datetime
import os
import numpy as np
import subprocess 

print('get it girl, and by girl i mean walrus')
resdir = '/data/tjarniko/results/may10_a7'
subprocess.call('ncks -v dissolved_inorganic_carbon SKOG_1h_20160331_20160414_ptrc_T.nc DIC_1h.nc', shell=True,cwd=resdir)
print('a7')
subprocess.call('ncks -v total_alkalinity SKOG_1h_20160331_20160414_ptrc_T.nc TA_1h.nc', shell=True,cwd=resdir)
print('a7')
subprocess.call('ncks -v dissolved_oxygen SKOG_1h_20160331_20160414_ptrc_T.nc OXY_1h.nc', shell=True,cwd=resdir)
print('a7')
subprocess.call('ncks -v vosaline SKOG_1h_20160331_20160414_grid_T.nc sn_1h.nc', shell=True,cwd=resdir)
print('a7')
subprocess.call('ncks -v votemper SKOG_1h_20160331_20160414_grid_T.nc tn_1h.nc', shell=True,cwd=resdir)
print('a7')

resdir = '/data/tjarniko/results/may10_a8'
subprocess.call('ncks -v dissolved_inorganic_carbon SKOG_1h_20160415_20160429_ptrc_T.nc DIC_1h.nc', shell=True,cwd=resdir)
print('a7')
subprocess.call('ncks -v total_alkalinity SKOG_1h_20160415_20160429_ptrc_T.nc TA_1h.nc', shell=True,cwd=resdir)
print('a7')
subprocess.call('ncks -v dissolved_oxygen SKOG_1h_20160415_20160429_ptrc_T.nc OXY_1h.nc', shell=True,cwd=resdir)
print('a7')
subprocess.call('ncks -v vosaline SKOG_1h_20160415_20160429_grid_T.nc sn_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v votemper SKOG_1h_20160415_20160429_grid_T.nc tn_1h.nc', shell=True,cwd=resdir)
print('a8')
