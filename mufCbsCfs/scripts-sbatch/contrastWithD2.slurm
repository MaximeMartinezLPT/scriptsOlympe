#!/bin/bash
#SBATCH -J chbd
#SBATCH -N 5
#SBATCH -n 180
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=05:00:00
#SBATCH --job-name="D2contrast"


module purge
module load openmpi/gnu/2.0.2.10 chdb/1.0-ompi python/3.6.3
ulimit -s 10240

nruns=180
name="Rectangle"
wdir=/tmpdir/p0110mm/mufCBSCFS/contrastWithD2-$name-${SLURM_JOB_ID}
inputfile=inputs/$name

python3 contrastWithD2.py initialize $wdir/ $inputfile $nruns
srun chdb --in-type "1 $nruns" --report "$wdir/report.txt" --command "python3 contrastWithD2.py compute $wdir/ %path%"
python3 contrastWithD2.py gather $wdir/
python3 contrastWithD2.py plot $wdir/
