#!/bin/bash
#SBATCH -J chbd
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=00:05:00
#SBATCH --job-name="qe(beta)"


module purge
module load openmpi/gnu/2.0.2.10 chdb/1.0-ompi python/3.6.3
ulimit -s 10240

nruns=35
name="DoublePuits-g0d25-e0d1"
wdir=/tmpdir/p0110mm/classical/PortraitDePhase-$name-${SLURM_JOB_ID}
inputfile=inputs/$name

python3 phasePortrait.py initialize $wdir/ $inputfile $nruns
srun chdb --in-type "1 $nruns" --report "$wdir/report.txt" --command "python3 phasePortrait.py compute $wdir/ %path%"
python3 phasePortrait.py gather $wdir/
python3 phasePortrait.py plot $wdir/
python3 phasePortrait.py plothusimi $wdir/
