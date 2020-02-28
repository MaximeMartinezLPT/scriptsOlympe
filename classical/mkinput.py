import sys
sys.path.insert(0, '..')
from inputparams import writeInput

dfile="1resclass"

params={}

params['tmax']=10000

params['gamma']=0.40
params['epsilon']=0.15
params['phi']=0.25

writeInput("inputs/"+dfile+".txt",params)


