import DAY_3panel_map as D3m


dirstr = './PRES_PLOT/'

#yvel
pdat = 'dic_ra'
print(pdat)
deep = 22
resdir = 'may10_a7/'
segment = 7
D3m.MRAP(resdir,dirstr,segment,pdat,deep)
resdir = 'may10_a8/'
segment = 8
D3m.MRAP(resdir,dirstr,segment,pdat,deep)
