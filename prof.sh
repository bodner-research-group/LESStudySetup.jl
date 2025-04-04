#!/bin/bash
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH -q regular
#SBATCH --account=m4672
#SBATCH --time=24:00:00
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

$JULIA --project --check-bounds=no -e 'using Pkg; Pkg.instantiate()'
srun ./launch.sh $JULIA --check-bounds=no --project experiments/benchmark_nonhydrostatic.jl 
