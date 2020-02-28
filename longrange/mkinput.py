import sys
sys.path.insert(0, '..')
from inputparams import writeInput

dfile="1resonance"

params={}

params['h']=0.40
params['gamma']=0.355
params['epsilon']=0.15
params['phi']=0.0

# ~ params['Ncell']=125
params['Npcell']=32




writeInput("inputs/"+dfile+".txt",params)


