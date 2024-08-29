module load cray-mpich

export SBATCH_ACCOUNT=m4499
export SALLOC_ACCOUNT=m4499

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_CUDA_MEMORY_POOL=none
export JULIA_NVTX_CALLBACKS=gc

export SLURM_CPU_BIND="cores"
export CRAY_ACCEL_TARGET="nvidia80"

module load julia

export JULIA="julia"
