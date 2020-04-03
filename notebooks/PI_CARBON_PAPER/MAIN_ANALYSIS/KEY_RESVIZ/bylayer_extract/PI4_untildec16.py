
#PI4
import importlib
import extract_bylayer_mean as ebm
importlib.reload(ebm)




#constant
start = '2015-01-01' #start of timeperiod
end = '2015-12-16' #end of timeperiod (typically a year)
inletmask = False #are we masking out Toba/Bute/Jervis?
sdir = '/PILA4/PI4' #where under directory tree do we find ncs 

#changes
shortdesc = 'PI4_DIC'
ftype = 'carp' #type of model result .nc 
varname = 'dissolved_inorganic_carbon' #name of variable
fname = shortdesc + '_means_inletsIN' #name of resulting pkl 
ebm.extractor(start, end, ftype, sdir, varname, fname,  inletmask)

shortdesc = 'PI4_TA'
ftype = 'carp' #type of model result .nc 
varname = 'total_alkalinity' #name of variable
fname = shortdesc + '_means_inletsIN' #name of resulting pkl 
ebm.extractor(start, end, ftype, sdir, varname, fname,  inletmask)

shortdesc = 'PI4_sal'
ftype = 'grid_T' #type of model result .nc 
varname = 'vosaline' #name of variable
fname = shortdesc + '_means_inletsIN' #name of resulting pkl 
ebm.extractor(start, end, ftype, sdir, varname, fname,  inletmask)

shortdesc = 'PI4_temp'
ftype = 'grid_T' #type of model result .nc 
varname = 'votemper' #name of variable
fname = shortdesc + '_means_inletsIN' #name of resulting pkl 
ebm.extractor(start, end, ftype, sdir, varname, fname,  inletmask)

shortdesc = 'PI4_nit'
ftype = 'ptrc' #type of model result .nc 
varname = 'nitrate' #name of variable
fname = shortdesc + '_means_inletsIN' #name of resulting pkl 
ebm.extractor(start, end, ftype, sdir, varname, fname,  inletmask)

shortdesc = 'PI4_diat'
ftype = 'ptrc' #type of model result .nc 
varname = 'diatoms' #name of variable
fname = shortdesc + '_means_inletsIN' #name of resulting pkl 
ebm.extractor(start, end, ftype, sdir, varname, fname,  inletmask)

