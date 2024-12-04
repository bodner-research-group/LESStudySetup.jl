#!/bin/bash
#SBATCH -C gpu
#SBATCH -N 256 
#SBATCH -q regular
#SBATCH --account=m4499
#SBATCH --time=10:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=none

source setup_perlmutter.sh

export SLURM_CPU_BIND="cores"
export CRAY_ACCEL_TARGET="nvidia80"

cat > launch.sh << EoF_s
#! /bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3
exec \$*
EoF_s
chmod +x launch.sh

$JULIA --check-bounds=no --project -e 'using Pkg; Pkg.instantiate()'
srun ./launch.sh $JULIA --check-bounds=no --project experiments/nonhydrostatic_experiment.jl 
