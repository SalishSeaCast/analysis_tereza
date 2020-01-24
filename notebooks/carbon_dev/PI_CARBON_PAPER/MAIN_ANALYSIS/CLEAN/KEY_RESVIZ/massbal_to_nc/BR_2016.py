from __future__ import print_function
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import scipy as sp
import seawater
import datetime as dt
""
from salishsea_tools import (
    nc_tools,
    viz_tools,
    geo_tools,
    tidetools
)
import netCDF4 as nc
import cmocean as cm
import glob
import sys
sys.path.append('/data/tjarniko/mocsy')
sys.path.append('/data/tjarniko/MEOPAR/at3/notebooks/carbon_dev/CCCmaDEV/CCCma_src')
import mocsy
import CCCma
import CCCma_stations as cs
from matplotlib import reload
import arrow
import gsw
import MASSBAL_ncmaker as mb

ncname_BR = 'MASSBAL_BR_2016_1_fullyear.nc'
sdir = 'BR_2016'
start = '2016-01-01'
end = '2016-12-31'
st = dt.datetime(2016,1,1)
en = dt.datetime(2016,12,31)
calc_JS = 1
mb.create_massbal_nc(ncname_BR, sdir, start, end, st, en, calc_JS)