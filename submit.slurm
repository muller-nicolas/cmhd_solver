#!/bin/sh
#SBATCH --account=c1715642
#SBATCH --job-name=test
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=00:10:00

##SBATCH --dependency=singleton
##SBATCH --mail-user=nicolas.muller@phys.ens.fr --mail-type=END

module purge
module load craype-x86-trento
module load craype-network-ofi
module load craype/2.7.30
module load intel-oneapi/2022.1.0
module load cray-mpich/8.1.26
module load cray-fftw/3.3.10.6

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lus/home/CT5/c1715642/muller/usr/lib
#ulimit -s unlimited
#OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

date
srun ./CMHD2D-v8
srun ./spectrum-anim
srun ./spectrum-anim2
srun ./spectrum-anim3
date
