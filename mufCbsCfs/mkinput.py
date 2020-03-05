import sys
import numpy as np

sys.path.insert(0, '..')
from inputparams import writeInput

dfile="Rectangle"

params={}
# ~ params['nDq']=25

# 1.Nombre de valeurs du paramètre varié pour calculer les D2
params['potential']="Rectangle"

# 2. Pour Square, largeur du potentiel
# ~ params['alpha']=0.5

# 3. Pour propagation temporelle
# ~ params['N']=1024
# ~ params['naverage']=32
# ~ params['i0']=int(1024/4)


writeInput("inputs/"+dfile+".txt",params)


