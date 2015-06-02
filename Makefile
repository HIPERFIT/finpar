# Supposed to be the Centralised makefile for benchmarks
# BUT IT IS NOT IMPLEMENTED RIGHT NOW. RATHER GO AND MAKE
# EACH BENCHMARK (VERSION) INDIVIDUALLY IN ITS OWN FOLDER,
# PLEASE READ the README.md file!!!!!
#
# This makefile is common infrastructure to compile and run all benchmarks.
#
# It recursively calls make in specified subdirectories with the
# benchmark programs, and expects that each such directory contains a
# makefile with the following targets:
#
# <default> compile the benchmark 
# run_X     run the benchmark with data set X
#            (where X is in {small,medium,large})
# clean
#
# Some centralised targets are defined here:
# run_all: runs all benchmarks with all data
# clean: calls make clean in all configured subdirectories
# OptionPricing: compile all programs in OptionPricing subdir
# 


#################################################################

# As per today, we have the following working directories: 

# Option Pricing
# ./OptionPricing/OrigCpp      -- sequential, original C(++) code
# ./OptionPricing/CppOpenMP    -- OpenMP version
# ./OptionPricing/CppOpenCL    -- GPU version using OpenCL
# ./OptionPricing/HaskellLH    -- Haskell version documenting all parallelism.
#
# Local Volatility Calibration
# ./LocVolCalib/OrigCpp        -- sequential, original C(++) code
# ./LocVolCalib/COpenMP        -- OpenMP C version of the code
# ./LocVolCalib/AllParOpenCLMP -- parallelizing the outer two loops in OpenCL and OpenMP
# ./LocVolCalib/OutParOpenCLMP -- parallelizing the whole loop nest in OpenCL and OpenMP
# ./LocVolCalib/HaskellLH      -- Haskell version documenting all parallelism.
#
# Interest Rate Calibration:
# ./InterestCalib/OrigCpp      -- sequential C++ code (translated from the original Caml code)
# ./InterestCalib/CppOpenMP    -- OpenMP version, in which only the outermost loop is parallel
# ./InterestCalib/CppOpenCL    -- OpenCL version, which needs to exploit parallelism on all levels.
# ./InterestCalib/HaskellLH    -- Haskell version documenting all parallelism.
#

BENCHMARKS =OptionPricing #LocVolCalib #InterestCalib

VERSIONS_OptionPricing = OptionPricing/OrigCpp OptionPricing/CppOpenMP \
			OptionPricing/HaskellLH OptionPricing/CppOpenCL

VERSIONS_LocVolCalib   = LocVolCalib/OrigCpp LocVolCalib/COpenMP LocVolCalib/AllParOpenCLMP \
			LocVolCalib/HaskellLH LocVolCalib/OutParOpenCLMP

VERSIONS_InterestCalib = InterestCalib/OrigCpp  InterestCalib/CppOpenMP \
			InterestCalib/CppOpenCL InterestCalib/HaskellLH

include platform.mk

####################### rules start here ########################

# default: help text
.PHONY: help
help:
	@echo "\tFor now, please read the makefile before using it"
	@echo "\t XXX Help text should be added soon here"

all	: $(BENCHMARKS) run_all

run_all ::
	@echo "Running all benchmarks"

OptionPricing: $(VERSIONS_OptionPricing)

LocVolCalib: $(VERSIONS_LocVolCalib)

InterestCalib: $(VERSIONS_InterestCalib)


################## target construction functions ################
#<name>	:
#	$(MAKE) -C <name>

#run_<name>	: <name>
#	make -C <name> run_<X>

# XXX BTW some renaming/reducing dir.s will be appropriate...
#  */Data/[small|medium|large][in|out] is good enough

define runRule
# arg.s: benchmark/version, Category (small,medium,large)
#$$(warning runRule($(1),$(2)))

.PHONY:	run_$(1)_$(2)
run_$(1)_$(2): $(1)
	$$(MAKE) -C $(1) run_$(2)
run_all :: run_$(1)_$(2)

endef

define mkStdRule
# arg.: Benchmark/Version#
#$$(warning mkStdRule($(1)))

.PHONY: $(1)
# building the target
$(1):
	$$(MAKE) -C $(1)

# running small,medium,large
$$(eval $$(call runRule,$(1),small))
$$(eval $$(call runRule,$(1),medium))
$$(eval $$(call runRule,$(1),large))

# and cleaning...
clean ::
	$$(MAKE) -C $(1) clean
endef

################## constructing the targets ################

$(foreach B,$(VERSIONS_OptionPricing),$(eval $(call mkStdRule,$(B))))
$(foreach B,$(VERSIONS_LocVolCalib),$(eval $(call mkStdRule,$(B))))

$(foreach B,$(VERSIONS_InterestCalib),$(eval $(call mkStdRule,$(B))))

# custom rules would go here as well
