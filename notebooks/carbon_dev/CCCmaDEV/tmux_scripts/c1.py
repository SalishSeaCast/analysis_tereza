import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import matplotlib.patches as patches
plt.style.use('seaborn-whitegrid')
import netCDF4 as nc
import cmocean as cm
import numpy as np
from salishsea_tools import (
    viz_tools,
)
import sys
sys.path.append('/data/tjarniko/mocsy')
import mocsy
import CCCma
import CCCma_stations as cs
#from matplotlib import reload
import arrow

start = '2016-01-01'
end = '2016-12-31'
CCCma.range_analyzer(start,end, prof = False, ncb = False, surfmap = False, plume = False, pspace = True)

