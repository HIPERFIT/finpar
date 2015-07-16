################################################
### Overview of the Workings of finpar ###
################################################

finpar contains three problems: interest rate calibration, local volume calibration, and option pricing. They each have a number of implementations. Each implementation is built manually inside the folder of /<problem>/<implementation>. Many benchmarks utilize OpenCL. The OpenCL implementation is (to some degree) configured in /platform.mk -- e.g., the number of cores in the graphic card and the number of cores in the cpu is set here.

Libraries which can be utilized for all implementations are located in /include/. Libraries which are only used by one implementation are located in <problem>/includeC.

The library functions generally contain a lot of preprocessor code and the function declerations are made in files with .h extensions. All library functions are defined in files with a .h extension. Whether this is standard procedure or not, I do not know.

Each problem contains a data folder at /<problem>/Data. For each problem, this folder contains three different pairs of inputs and outputs: a computationally easy pair, a medium pair, and a large pair with the respective names: Small, Medium, Large. When the simulation is run, the algorithm initializes the data with that of the input file. The calculated output is then compared with the result stored in the output file and the test validates if the results are the same (or within some margin of error since floats are used?). Upon completion, the test also shows the elapsed time.

So each problem can have several different implementations which utilize the same input/output data. But the make files are individual for each implementation. It should not be necessary to have seperate make files since the compiling and running of the programs should be managed centrally.

################################################
### Changes Made to This Branch ###
################################################

The Makefiles for each implementatation and the SDK_stub.h have been modified such that the compilation (building of kernel) receives the argument "-I <path of invoked Makefile>". This occurs in this function call:
ciErr1 = clBuildProgram(cpProgram, 1, cdDevices+dev_id, newCompileOptions, NULL, NULL);

where "newCompileOptions" includes the above mentioned argument. This makes it possible to run on non-nVidia platforms since the compilation may occur somewhere else than the location of the Makefile, making the compiler unable to locate any resources that needs to be utilized, had this fix not been applied.

