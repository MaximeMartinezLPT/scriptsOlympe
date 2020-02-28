#!/bin/bash
#SBATCH -J chbd
#SBATCH -N 5
#SBATCH -n 180
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=00:02:00
#SBATCH --job-name="qe(beta)"


module purge
module load openmpi/gnu/2.0.2.10 chdb/1.0-ompi python/3.6.3
ulimit -s 10240

nruns=179
name="1resonance"
wdir=/tmpdir/p0110mm/new/longrange/spectrum-beta-$name-${SLURM_JOB_ID}
inputfile=inputs/$name

python3 spectrum-beta.py initialize $wdir/ $inputfile $nruns
srun chdb --in-type "1 $nruns" --report "$wdir/report.txt" --command "python3 spectrum-beta.py compute $wdir/ %path%"
python3 spectrum-beta.py gather $wdir/
python3 spectrum-beta.py plot $wdir/
