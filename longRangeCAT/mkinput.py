import sys
sys.path.insert(0, '..')
from inputparams import writeInput

dfile="g0d20-e0d15-h0d45"

params={}

params['h']=0.45
params['gamma']=0.20
params['epsilon']=0.15
params['phi']=0.0

# ~ params['Ncell']=125
params['Npcell']=32




writeInput("inputs/"+dfile+".txt",params)


