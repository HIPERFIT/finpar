# Centralised makefile for benchmarks
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
# GenericPricing: compile all programs in GenericPricing subdir
# 


#################################################################

# As per today, we have the following working directories: 

# ./GenericPricing/Orig_COpenMp
# ./GenericPricing/HaskellLH
# ./GenericPricing/CppOpenCL
#
# planned future work:
# ./CalibVolDiff/Orig_COpenMP
# ./CalibVolDiff/VectAll
# ./CalibVolDiff/Original
# ./CalibVolDiff/VectOuters
# ./CalibGA/CppAndGPU
# ./CalibGA/python
# ./CalibGA/OCaml

BENCHMARKS =GenericPricing #CalibVolDiff #CalibGA

VERSIONS_GenericPricing =GenericPricing/Orig_COpenMp GenericPricing/HaskellLH GenericPricing/CppOpenCL

VERSIONS_CalibVolDiff   =CalibVolDiff/Orig_COpenMP CalibVolDiff/VectAll \
                         CalibVolDiff/Original CalibVolDiff/VectOuters
VERSIONS_CalibGA        =CalibGA/CppAndGPU CalibGA/python CalibGA/OCaml

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

GenericPricing: $(VERSIONS_GenericPricing)

#CalibVolDiff: $(VERSIONS_CalibVolDiff)
#CalibGA: $(VERSIONS_CalibGA)


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

$(foreach B,$(VERSIONS_GenericPricing),$(eval $(call mkStdRule,$(B))))
# $(foreach B,$(VERSIONS_CalibGA),$(eval $(call mkStdRule,$(B))))
# $(foreach B,$(VERSIONS_CalibVolDiff),$(eval $(call mkStdRule,$(B))))

# custom rules would go here as well
