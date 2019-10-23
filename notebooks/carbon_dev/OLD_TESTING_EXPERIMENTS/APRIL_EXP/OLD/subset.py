import netCDF4 as nc
import datetime
import os
import numpy as np
import subprocess 

print('starting subset')
resdir = '/data/tjarniko/results/may10_a1'
subprocess.call('ncks -v votemper SKOG_1h_20160101_20160115_grid_T.nc tn_1h.nc', shell=True,cwd=resdir)
print('walrus a1')

resdir = '/data/tjarniko/results/may10_a2'
subprocess.call('ncks -v dissolved_inorganic_carbon SKOG_1h_20160116_20160130_ptrc_T.nc DIC_1h.nc', shell=True,cwd=resdir)
print('walrus')
subprocess.call('ncks -v total_alkalinity SKOG_1h_20160116_20160130_ptrc_T.nc TA_1h.nc', shell=True,cwd=resdir)
print('walrus')
subprocess.call('ncks -v dissolved_oxygen SKOG_1h_20160116_20160130_ptrc_T.nc OXY_1h.nc', shell=True,cwd=resdir)
print('walrus')
subprocess.call('ncks -v vosaline SKOG_1h_20160116_20160130_grid_T.nc sn_1h.nc', shell=True,cwd=resdir)
print('walrus')
subprocess.call('ncks -v votemper SKOG_1h_20160116_20160130_grid_T.nc tn_1h.nc', shell=True,cwd=resdir)
print('a2')

resdir = '/data/tjarniko/results/may10_a3'
subprocess.call('ncks -v dissolved_inorganic_carbon SKOG_1h_20160131_20160214_ptrc_T.nc DIC_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v total_alkalinity SKOG_1h_20160131_20160214_ptrc_T.nc TA_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v dissolved_oxygen SKOG_1h_20160131_20160214_ptrc_T.nc OXY_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v vosaline SKOG_1h_20160131_20160214_grid_T.nc sn_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v votemper SKOG_1h_20160131_20160214_grid_T.nc tn_1h.nc', shell=True,cwd=resdir)
print('a3')

resdir = '/data/tjarniko/results/may10_a4'
subprocess.call('ncks -v dissolved_inorganic_carbon SKOG_1h_20160214_20160229_ptrc_T.nc DIC_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v total_alkalinity SKOG_1h_20160214_20160229_ptrc_T.nc TA_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v dissolved_oxygen SKOG_1h_20160214_20160229_ptrc_T.nc OXY_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v vosaline SKOG_1h_20160214_20160229_grid_T.nc sn_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v votemper SKOG_1h_20160214_20160229_grid_T.nc tn_1h.nc', shell=True,cwd=resdir)
print('a4')

resdir = '/data/tjarniko/results/may10_a5'
subprocess.call('ncks -v dissolved_inorganic_carbon SKOG_1h_20160301_20160315_ptrc_T.nc DIC_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v total_alkalinity SKOG_1h_20160301_20160315_ptrc_T.nc TA_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v dissolved_oxygen SKOG_1h_20160301_20160315_ptrc_T.nc OXY_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v vosaline SKOG_1h_20160301_20160315_grid_T.nc sn_1h.nc', shell=True,cwd=resdir)
subprocess.call('ncks -v votemper SKOG_1h_20160301_20160315_grid_T.nc tn_1h.nc', shell=True,cwd=resdir)
print('a5')
