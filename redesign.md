(To eliminate ambiguity, here is the nomenclature: we have a number of
*benchmarks* (currently CalibGA, CalibVolDiff, and GenericPricer),
each of which have several *data sets* (typically Small, Medium, and
Large), and several *implementations* (right now mostly different
versions of C) each of which may have several *configurations*.
Running a benchmark consists of selecting a data set and an
implementation, and possibly specifying a specific configuration of
the implementation.

Recently, a Martin, Frederik, and myself have been implementing the
finpar bechmarks in more diverse programming languages - (streaming)
NESL, APL and Futhark, at least.  Unfortunately, the current finpar
infrastructure is not very easy to work with, and so their work has
not been integrated.  I have identified the following problems:

  * Implementations are not cleanly separated from data sets and
    ancillary code.

    **Solution**: for each benchmark, have a directory that contains
    only implementations.

  * Building an implementation modifies the implementation
    directory, and more importantly, configuring an implementation
    often involves manually modifying files in the directory (see
    CalibGA/includeC/KerConsts.h for an example).  This is really
    bad and makes structured and reproducible benchmarking almost
    impossible.

    **Solution**: when "compiling" an implementation, put everything
    in a new, separate directory, I will call the *instantiation
    directory*.  All configuration must be done via passing options to
    the compilation step, and will be reflected in the files put in
    the instantiation directory.

  * Adding new implementations is a mess, because you have to modify
    the global build system.

    **Solution**: define a setup/run-protocol that each benchmark
    implementation must follow, and which can be used by a generic
    controller script.

  * Validation is done by the benchmark implementations.  There is no
    reason to do this.

    **Solution**: have the implementation produce
    their results in some well-defined format, and have the controller
    script validate it.

  * Everything is done with Makefiles.  Nobody likes modifying
    Makefiles, and we don't need incremental rebuilds anyway.


    **Solution**: write as much as possible in Python or simple shell
    script.

I propose the following rough protocol:

  * Each benchmark implementation must include one executable file,
    called `instantiate`.  This can be written in whatever language
    one prefers.

  * When the `instantiate` program for an implementation is invoked,
    the following environment variables must be set:

    * `FINPAR_IMPLEMENTATION`, which must point at the
      implementation directory.  This is to get around the fact that
      it's not always easy to find the location of the running
      program.

    * `FINPAR_DATASET`, which must point at a directory containint
    `.input` and `.output` files.

  * The `instantiate` program will instantiate the implementation in
    the *current directory*, which will become the instantiation
    directory.

  * The `instantiate` program can be passed command-line options to
    futher configure the implementation.  These are defined on a
    per-implementation basis, and not standardised.

  * After instantiation, the instantiation directory must contain a
    program `run`, which, when executed, will run the benchmark
    implementation.  The result will be two files in the
    instantiation directory:

    * `runtime.txt`, which contains the runtime in milliseconds as
      an integer.

    * `result.data`, which contains the result in ou well-defined
      data format.

    I have decided that the runtime should be measured by the
    implementation itself, as it is not possible to black-box
    measure this without possibly measuring the wrong things (like
    kernel compilation, exotic hardware setup, parsing of input
    data, or IO).

The following questions have yet to be answered:

  * What data format should we use?  Currently, finpar uses the
    Futhark value format, which is pretty simple.  We can possibly
    make life even simpler by using JSON, but it is incredibly
    annoying that JSON does not support comments.

  * Should we use environment variables at all?  It was mostly to
    avoid having the instantiate-script do command-line parsing unless
    it wants to.

Yet, I think this is a good protocol.  It will allow us to build an
easy-to-use controller program on top of it, that can automatically
generate a bunch of different instantiations with different
configurations and data sets, and maybe draw graphs of the results,
etc.  I estimate that the above could be implemented fairly quickly,
and sanity-checked with the extant benchmark implementations.
