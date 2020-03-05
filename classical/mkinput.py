import sys
sys.path.insert(0, '..')
from inputparams import writeInput

dfile="DoublePuits-g0d25-e0d1"

params={}

# ~ params['tmax']=10000
params['tmax']=5000

params['gamma']=0.25
params['epsilon']=1
params['phi']=0.00

writeInput("inputs/"+dfile+".txt",params)


