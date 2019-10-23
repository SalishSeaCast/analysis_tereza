import phys_map as pm

dirstr = './APR_velar/'

#
segment = 4
resdir = 'may10_a4/'
zlevel = 0
figtit = 'VA_zl_' + str(zlevel) + '_hr_'
pm.MRAP(resdir,dirstr,figtit,segment,zlevel)
# 34
zlevel = 22
figtit = 'VA_zl_' + str(zlevel) + '_hr_'
pm.MRAP(resdir,dirstr,figtit,segment,zlevel)
# 75
zlevel = 25
figtit = 'VA_zl_' + str(zlevel) + '_hr_'
pm.MRAP(resdir,dirstr,figtit,segment,zlevel)

#
segment = 5
resdir = 'may10_a5/'
zlevel = 0
figtit = 'VA_zl_' + str(zlevel) + '_hr_'
pm.MRAP(resdir,dirstr,figtit,segment,zlevel)
# 34
zlevel = 22
figtit = 'VA_zl_' + str(zlevel) + '_hr_'
pm.MRAP(resdir,dirstr,figtit,segment,zlevel)
# 75
zlevel = 25
figtit = 'VA_zl_' + str(zlevel) + '_hr_'
pm.MRAP(resdir,dirstr,figtit,segment,zlevel)
