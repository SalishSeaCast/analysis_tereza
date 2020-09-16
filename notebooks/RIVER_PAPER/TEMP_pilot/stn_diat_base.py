import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
plt.style.use('classic')
from matplotlib.patches import Rectangle
import pickle
from cmocean import cm
import warnings
import netCDF4 as nc
import arrow
import glob
import sys
sys.path.append('../RIVERS_pilot/bylayer_extract/')
import extract_bylayer_mean_BOXMASK as ebmBOX 
import extract_bylayer_mean as ebm
from importlib import reload
from salishsea_tools import viz_tools


warnings.simplefilter('ignore')

plt.rcParams.update({'font.size': 13,
                     'xtick.labelsize' : 13,
                     'ytick.labelsize' : 13})
                     
reload(ebmBOX)

reload(ebmBOX)

jdf_x = 85; jdf_y = 325; jdf_col = 'teal'
jdf2_x = 50; jdf2_y = 370; jdf2_col = 'mediumaquamarine'
jdf3_x = 130; jdf3_y = 290; jdf3_col = 'lightseagreen'

cs_x = 250; cs_y = 500; cs_col = 'royalblue'
cs2_x = 280; cs2_y = 440; cs2_col = 'navy'
cs3_x = 220; cs3_y = 540; cs3_col = 'dodgerblue'

haro_x = 230; haro_y = 310; haro_col = 'tomato'
haro2_x = 340; haro2_y = 290; haro2_col = 'firebrick'
haro3_x = 260; haro3_y = 335; haro3_col = 'indianred'


ns_x = 160; ns_y = 680; ns_col = 'olive'
ns2_x = 150; ns2_y = 640; ns2_col = 'yellowgreen'
ns3_x = 155; ns3_y = 710; ns3_col = 'palegoldenrod'

#(start, end, ftype, sdir, varname, fname, y, x)
pkldir = './pkls/'
start = '2017-01-01' #start of timeperiod
end = '2017-12-31' #end of timeperiod (typically a year)
ftype = 'ptrc' #type of model result .nc 
sdir = 'BASE/ncs/' #where under directory tree do we find ncs 
varname = 'diatoms' #name of variable


fname = pkldir + 'BASE_diat_means_cs2_BOX'
y = cs2_y; x = cs2_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_cs_BOX'
y = cs_y; x = cs_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_cs3_BOX'
y = cs3_y; x = cs3_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_haro2_BOX'
y = haro2_y; x = haro2_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_haro_BOX'
y = haro_y; x = haro_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_haro3_BOX'
y = haro3_y; x = haro3_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_jdf2_BOX'
y = jdf2_y; x = jdf2_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_jdf_BOX'
y = jdf_y; x = jdf_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_jdf3_BOX'
y = jdf3_y; x = jdf3_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_ns2_BOX'
y = ns2_y; x = ns2_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_ns_BOX'
y = ns_y; x = ns_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
fname = pkldir + 'BASE_diat_means_ns3_BOX'
y = ns3_y; x = ns3_x
ebmBOX.extractor(start, end, ftype, sdir, varname, fname, y, x )
