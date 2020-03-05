#!/bin/bash
#SBATCH -J chbd
#SBATCH -N 5
#SBATCH -n 180
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=00:02:00
#SBATCH --job-name="ebeta"


module purge
module load openmpi/gnu/2.0.2.10 chdb/1.0-ompi python/3.6.3
ulimit -s 10240

nruns=179
name="g0d20-e0d15-h0d45"
wdir=/tmpdir/p0110mm/longRangeCAT/spectrumBeta-$name-${SLURM_JOB_ID}
inputfile=inputs/$name

python3 qenergiesWithBeta.py initialize $wdir/ $inputfile $nruns
srun chdb --in-type "1 $nruns" --report "$wdir/report.txt" --command "python3 qenergiesWithBeta.py compute $wdir/ %path%"
python3 qenergiesWithBeta.py gather $wdir/
python3 qenergiesWithBeta.py plot $wdir/
