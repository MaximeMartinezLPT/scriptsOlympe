#!/bin/bash
#SBATCH -J chbd
#SBATCH -N 3
#SBATCH -n 108
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=00:25:00
#SBATCH --job-name="qe(beta)"


module purge
module load openmpi/gnu/2.0.2.10 chdb/1.0-ompi python/3.6.3
ulimit -s 10240

nruns=107
name=g0d40-e0d50
wdir=/tmpdir/p0110mm/new/longrange/spp-$name-${SLURM_JOB_ID}
inputfile=inputs/$name

python3 stroboscopicPhasePortrait.py initialize $wdir/ $inputfile $nruns
srun chdb --in-type "1 $nruns" --report "$wdir/report.txt" --command "python3 stroboscopicPhasePortrait.py compute $wdir/ %path%"
python3 stroboscopicPhasePortrait.py gather $wdir/
python3 stroboscopicPhasePortrait.py plot $wdir/
