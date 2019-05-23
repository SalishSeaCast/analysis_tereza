import glob
import subprocess

for fn in glob.iglob('/data/tjarniko/results/BR_2nd_2015_cop/SKOG_2/ncs/*.nc'):
    print(fn)
    output_fn = fn
    
    subprocess.run(['ncks', '-4', '-L4', '-O', fn, output_fn])
