# platform specification for benchmarks. Edit as appropriate for your
# platform before running

# 1 = platform has a GPU
HAVE_GPU      = 1
GPU_DEVICE_ID = 0

# GPU warp size. Keep consistent if you use the non-log variant
#GPU_WARP      = 32
GPU_LG_WARP   = 4

# GPU memory sizes in kilobyte
GPU_LOCAL_MEM = 48
GPU_CONST_MEM = 64
GPU_REG_MEM   = 64
# device memory in gigabyte
GPU_DEVICE_MEM= 1
# ``Optimal'' Amount of Local/Fast Memory Per Thread 
GPU_LOCAL_MEM_PER_TH=8
# Number of GPU cores
GPU_NUM_CORES = 64

# CPU and memory spec.
NCORES = 32
# in gigabyte
MEMORY = 1

export
