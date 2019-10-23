import DAY_yvel_map as Dym

dirstr = './APR_yvel_day/'
figtit = 'YVEL_day_'



segment = 6
print('segment' + str(segment))
resdir = 'may10_a6/'
Dym.MRAP(resdir,dirstr,figtit,segment)

segment = 7
print('segment' + str(segment))
resdir = 'may10_a7/'
Dym.MRAP(resdir,dirstr,figtit,segment)

segment = 8
print('segment' + str(segment))
resdir = 'may10_a8/'
Dym.MRAP(resdir,dirstr,figtit,segment)