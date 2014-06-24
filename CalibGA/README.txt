1. Compile & Run with:

        $ make clean; make; make run_medium

2. One may change options 

        WITH_SOBOL, GPU_VERSION, and WITH_FLOAT

   in file KerConsts.h 

3. Benchmark Description:

A dynamic evolution model method, i.e., genetic algorithm,
for calibrating the interest rate based on a known history
of swaption prices.   

Briefly, the interest rate is modelled as a sum of two stochastic 
processes, which gives four unknown (real) parameters, and in 
addition the two processes are assumed correlated as well, 
i.e., a fifth parameter.

These five (unknown) parameters appear in the formula that 
computes the swaption's price, i.e., numerical integration
via hermitian-polynomials approximation.

The genetic algorithm is used to find the five parameters 
that best fit the (known) history of swaption prices.
 
