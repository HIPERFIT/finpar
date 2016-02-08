################################################
### Library of Financial Parallel Benchmarks ###
################################################

As per today, we have the following three benchmarks implemented in the
working directories: 

1. Option Pricing
  ./OptionPricing/OrigCpp      -- sequential, original C(++) code
  ./OptionPricing/CppOpenMP    -- OpenMP version
  ./OptionPricing/CppOpenCL    -- GPU version using OpenCL. Various optimizations, such as memory coalescing, fussion-fission, Sobol strength reduction, can be (de)selected by macro definition in file "./OptionPricing/includeC/Optimizations.h".
  ./OptionPricing/HaskellLH    -- Haskell version documenting all parallelism.

2. Local Volatility Calibration
  ./LocVolCalib/OrigCpp        -- sequential, original C(++) code
  ./LocVolCalib/COpenMP        -- OpenMP C version of the code
  ./LocVolCalib/AllParOpenCLMP -- parallelizing the outer two loops in OpenCL and OpenMP
  ./LocVolCalib/OutParOpenCLMP -- parallelizing the whole loop nest in OpenCL and OpenMP
  ./LocVolCalib/HaskellLH      -- Haskell version documenting all parallelism.

3. Interest Rate Calibration:
  ./InterestCalib/OrigCpp      -- sequential C++ code (translated from the original Caml code)
  ./InterestCalib/CppOpenMP    -- OpenMP version, in which only the outermost loop is parallel
  ./InterestCalib/CppOpenCL    -- OpenCL version, which needs to exploit parallelism on all levels.
  ./InterestCalib/HaskellLH    -- Haskell version documenting all parallelism.

#########################################################################
### Hardware/software description and Where was the benchmark tested: ###
#########################################################################

1. Please fill in various env variables in file "setup.mk".

2. Please describe the GPU and CPU hardware in file "platform.mk", e.g.,
    (i)  the size of GPU's global memory (GPU_DEVICE_MEM) or 
   (ii)  the number of cores that execute in lock step (SIMD), e.g,
         two to power GPU_LG_WARP (for NVIDIA, GPU_LG_WARP is 5 since 
         warp size is 32) or 
   (iii) the OpenCL identifier of the GPU device, e.g, GPU_DEVICE_ID = 0.

3. So far, the benchmark was tested mostly on NVIDIA GPUs.
   Please note that the current OpenCL implementation does not adhere to, 
   the OpenCL standard i.e., we have used (segmented) reduction and scan
   implementation that assumes the existance of SIMD instructions
   (of length "2 to power GPU_LG_WARP"), and intra-block barriers.
   (In our tests some INTEL/AMD GPUs required "GPU_LG_WARP == 4".)

   Specifically, the OpenCL implementation of "LocVolCalib/AllParOpenCLMP"
   uses intra-block segmented scan with the following operators:
    (i) 2x2 matrix multiplication, and
   (ii) linear-function composition.

4. It is assumed that gcc is installed (with support for OpenMP, i.e., -fopenmp)
    and an OpenCL implementation is available -- see file "setup.mk".

5. For the OpenMP version, one may set the desired number of OpenMP threads with:
        $ export OMP_NUM_THREADS = 8

#############################################
### How to Compile and Run the Benchmarks ###
#############################################

1. Currently, there is no over-arching script that allows to select, 
    make and run a specified set of benchmarks. Instead, one must go
    to one of the working directories and make and run the benchmark 
    there. For example, for the "Local Volatility Calibration" version
    that parallelizes all (parallel) loops, go to the corresponding folder
        $ cd LocVolCalib/AllParOpenCLMP

    then compile the GPU (OpenCL) or CPU (OpenMP) version:
        $ make clean; make gpu
    or
        $ make clean; make cpu

    then run it on the desired (small/medium/large) dataset:
        $ make run_small
        $ make run_medium
        $ make run_large

2. The input/reference datasets are located in folder "Data" of each benchmark
    (under the small, medium, and large subfolders).

3. For example to run the large dataset of the "Local Volatility Calibration"
    benchmark in which the two outermost loops have been parallelized, do:
        $ LocVolCalib/OutParOpenCLMP
        $ make clean; make gpu
        $ make run_large
    The output specifies 
    (i)   the version of the benchmark and dataset that is has been run
    (ii)  whether the result matches the reference result (1 for VALID and 0 for INVALID)
    (iii) the runtime in microseconds: 
            a) INCLUDES memory-input/output transfer and kernel-running time, but 
            b) EXCLUDES reading/wrting the dataset from/to a file, the validation 
               of the result and various (other) OpenCL queries, e.g., for the
               size of the GPU global/shared memory, etc.
    (iv)  the number of threads that have been used, and
     (v)  the array/scalar results. 
    For example the following is an output of the run:
        // GPU Outer-Level Massively-Parallel Execution of Volatility Calibration Benchmark:
        // OUTER=128, NUM_X=256, NUM_Y=256, NUM_T=64.
        1		// VALID   Result,
        857606  // Runtime in microseconds,
        1536	// GPU Threads,

         [ 0.029996 , 0.029204 , 0.028803 , 0.028405 , 0.028013 , 0.027626 , 0.027245 , 0.026870 , 0.026500 , 0.026136 , 0.025778 , 0.025426 , 0.025080 , 0.024739 , 0.024404 , 0.024075 , 0.023751 , 0.023433 , 0.023120 , 0.022813 , 0.022511 , 0.022214 , 0.021923 , 0.021636 , 0.021355 , 0.021078 , 0.020806 , 0.020538 , 0.020273 , 0.020014 , 0.019791 , 0.019512 , 0.019270 , 0.019033 , 0.018799 , 0.018568 , 0.018341 , 0.018117 , 0.017898 , 0.017682 , 0.017470 , 0.017261 , 0.017056 , 0.016855 , 0.016657 , 0.016463 , 0.016271 , 0.016083 , 0.015898 , 0.015716 , 0.015537 , 0.015361 , 0.015188 , 0.015018 , 0.014850 , 0.014686 , 0.014523 , 0.014364 , 0.014206 , 0.014052 , 0.013899 , 0.013749 , 0.013601 , 0.013456 , 0.013312 , 0.013171 , 0.013032 , 0.012895 , 0.012760 , 0.012626 , 0.012495 , 0.012365 , 0.012238 , 0.012112 , 0.011988 , 0.011865 , 0.011744 , 0.011625 , 0.011508 , 0.011392 , 0.011277 , 0.011164 , 0.011053 , 0.010942 , 0.010834 , 0.010726 , 0.010620 , 0.010515 , 0.010412 , 0.010309 , 0.010208 , 0.010108 , 0.010010 , 0.009912 , 0.009816 , 0.009720 , 0.009626 , 0.009533 , 0.009440 , 0.009349 , 0.009259 , 0.009170 , 0.009081 , 0.008994 , 0.008908 , 0.008822 , 0.008737 , 0.008653 , 0.008570 , 0.008488 , 0.008407 , 0.008326 , 0.008246 , 0.008167 , 0.008089 , 0.008011 , 0.007934 , 0.007858 , 0.007782 , 0.007707 , 0.007633 , 0.007560 , 0.007487 , 0.007414 , 0.007343 , 0.007271 , 0.007201 , 0.007131  ]	//Volatility Calibration Result


4. Note that most of the benchmarks versions support either the 
    OpenMP (CPU) or the OpenCL (GPU) implementations, i.e.,
    "make cpu" or "make gpu" -- see the working-directory description 
    above.  Some implementations support both CPU and GPU implementations,
    the rationale beeing that the provided CPU implementation is "similar"
    to the OpenCL one, i.e., the code has been transformed in the same
    way, e.g., array expansion, loop distribution and interchange.

    a) Each benchmark exhibits an "original", sequential version and
       an OpenMP version that are to be run on CPU. 

    b) The Haskell versions are aimed to fully specified application's
       parallelism. They use a very basic Haskell (map, reduce, filter, 
       scan on lists), so that they are "easy" to read. From this purpose
       they are also extremely slow: Why they (eventually) validate, they
       take a very long time to run. This is not necessarily because
       Haskell is slow; in particular efficient Haskell code can be
       written for every one of the benchmarks, but this has not been
       the goal of the current implementation.

    c) One or several OpenCL versions are of course provided for each bench.


5. Before compiling and running the benchmarks, check the sections 
    below to see if there is OS specific software that you need to 
    install before proceeding.

#########################
### MacOS Assumptions ###
#########################

* GCC 4.9 (non-clang version) is assumed as some of the benchmarks are
  making use of OpenMP. The Makefiles are assuming g++-4.9 to be
  accessible from the PATH.

    bash-3.2$ g++-4.9 --version
    g++-4.9 (GCC) 4.9.0 20140119 (experimental)
    Copyright (C) 2014 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

### Adding New Benchmarks

To add a new benchmark, first read the description of the execution
strategy in the Makefile located in the top-level directory.
