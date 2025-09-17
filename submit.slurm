#!/bin/sh
#SBATCH --account=c1715642
#SBATCH --job-name=nth1
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00

##SBATCH --dependency=singleton
##SBATCH --mail-user=nicolas.muller@lpp.polytechnique.fr --mail-type=END

module purge
module load PrgEnv-cray
module load craype-x86-genoa
module load cray-fftw

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

date
srun ./CMHD2D.exe
date
