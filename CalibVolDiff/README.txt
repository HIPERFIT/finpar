
1.  This is a benchmark received from Nordea, via Prof. Brian Vinter. 
    Here is Christian's (Andreetta) guess about what the application is about:
    " The task is stochastic volatility calibration, i.e., given a set of 
        (observed) prices of contracts, we identify the parameters of a model 
        of such prices, as a function of volatility (unknown), time and strikes (known), 
        and unobserved parameters like alpha, beta, nu, etc.

      In this case, the volatility is modelled as a system of continuous
        partial differential equations, which are solved via Crank-Nicolson's
        finite differences method.

      The model seems to be a variation of SABR.
    "

2. STRUCTURE:

    a) `NordeaOrig' folder contains the original sequential code, which can be simply
        compiled by: `$ g++ HiperfitExample1.cpp'

    b) `VectOuters' folder contains the GPU parallel version of the code, in which
        two (outer) dimensions of parallelism are exploited, but the inner dimension,
        corresponding to TRIDAG is sequentially executed. 

        The CPU execution executes in parallel the outermost loop (OpenMP).    
    
        The degree of parallelism is sensitive to the following parameters, set in
        `DataStructConst.h' file: OUTER_LOOP_COUNT, NUM_X, NUM_Y, which are all set
        by default to 128, but can be varied. The degree of exploited parallelism is:
            OUTER_LOOP_COUNT * MIN(NUM_X, NUM_Y).

        If you wish the result to not be printed, set DEBUG to 0 in DataStructConst.h,
        and if you wish GPU hardware info, define DEBUG_PRINT_GPU_INFO.

        To build for CPU execution, set OMP_NUM_THREADS, and:
            $make cpu   
        To build for GPU execution:
            $make gpu

        To run:
            $./nordea
        
    c) `VectAll' folder contains the GPU versions that exploits all parallel dimensions
        of parallelism, i.e., TRIDAG is parallelized as a combination of scan with
        2x2 matrix multiplication and linear function composition, respectively.

        It follows that it works well even with smaller values for the parallelism-degree
        sensitive parameters, e.g., OUTER_LOOP_COUNT x NUM_X x NUM_Y = 128 x 32 x 32.

        The rest is as with `VectOuters', e.g., in the CPU version the outermost loop
        (OUTER_LOOP_COUNT) is executed in parallel via OpenMP.

