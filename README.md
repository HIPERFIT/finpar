## Library of Financial Parallel Benchmarks

### How to Compile and Run the Benchmarks 

Before compiling and running the benchmarks, check the sections below
to see if there is OS specific software that you need to install
before proceeding.



### MacOS Assumptions

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
