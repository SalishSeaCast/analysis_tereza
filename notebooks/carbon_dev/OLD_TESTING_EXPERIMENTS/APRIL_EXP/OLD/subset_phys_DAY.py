import netCDF4 as nc
import datetime
import os
import numpy as np
import subprocess 

print('starting subset')
resdir = '/data/tjarniko/results/may10_a1'
subprocess.call('ncks -v vozocrtx SKOG_1d_20160101_20160115_grid_U.nc u_1d.nc', shell=True,cwd=resdir)
print('u')
subprocess.call('ncks -v time_counter SKOG_1d_20160101_20160115_grid_U.nc timecount_1d.nc', shell=True,cwd=resdir)
print('t')
subprocess.call('ncks -v vomecrty SKOG_1d_20160101_20160115_grid_V.nc v_1d.nc', shell=True,cwd=resdir)
print('v')

print('a2')

resdir = '/data/tjarniko/results/may10_a2'

subprocess.call('ncks -v vozocrtx SKOG_1d_20160116_20160130_grid_U.nc u_1d.nc', shell=True,cwd=resdir)
print('u')
subprocess.call('ncks -v time_counter SKOG_1d_20160116_20160130_grid_U.nc timecount_1d.nc', shell=True,cwd=resdir)
print('t')
subprocess.call('ncks -v vomecrty SKOG_1d_20160116_20160130_grid_V.nc v_1d.nc', shell=True,cwd=resdir)
print('v')

resdir = '/data/tjarniko/results/may10_a3'
subprocess.call('ncks -v vozocrtx SKOG_1d_20160131_20160214_grid_U.nc u_1d.nc', shell=True,cwd=resdir)
print('u')
subprocess.call('ncks -v time_counter SKOG_1d_20160131_20160214_grid_U.nc timecount_1d.nc', shell=True,cwd=resdir)
print('t')
subprocess.call('ncks -v vomecrty SKOG_1d_20160131_20160214_grid_V.nc v_1d.nc', shell=True,cwd=resdir)
print('v')

print('a3')

resdir = '/data/tjarniko/results/may10_a4'
subprocess.call('ncks -v vozocrtx SKOG_1d_20160214_20160229_grid_U.nc u_1d.nc', shell=True,cwd=resdir)
print('u')
subprocess.call('ncks -v time_counter SKOG_1d_20160214_20160229_grid_U.nc timecount_1d.nc', shell=True,cwd=resdir)
print('t')
subprocess.call('ncks -v vomecrty SKOG_1d_20160214_20160229_grid_V.nc v_1d.nc', shell=True,cwd=resdir)
print('v')

print('a4')

resdir = '/data/tjarniko/results/may10_a5'
subprocess.call('ncks -v vozocrtx SKOG_1d_20160301_20160315_grid_U.nc u_1d.nc', shell=True,cwd=resdir)
print('u')
subprocess.call('ncks -v time_counter SKOG_1d_20160301_20160315_grid_U.nc timecount_1d.nc', shell=True,cwd=resdir)
print('t')
subprocess.call('ncks -v vomecrty SKOG_1d_20160301_20160315_grid_V.nc v_1d.nc', shell=True,cwd=resdir)
print('v')

print('a5')
