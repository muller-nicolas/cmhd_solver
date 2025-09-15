#!/bin/sh
#SBATCH --account=c1715642
#SBATCH --job-name=test
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --time=01:00:00

##SBATCH --dependency=singleton
##SBATCH --mail-user=nicolas.muller@phys.ens.fr --mail-type=END

module purge
module load craype-x86-trento
module load craype-network-ofi
module load craype/2.7.30
module load intel-oneapi/2022.1.0
module load cray-mpich/8.1.26
module load cray-fftw/3.3.10.6

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

date
srun ./CMHD2D.exe
date
