FINPAR_BASE_DIR := $(shell pwd)
PROBLEM_NAMES := InterestCalib
INTERESTCALIB_DIRS := CppOpenMP OrigCpp
INTERESTCALIB_BIN_DIR := $(FINPAR_BASE_DIR)/bin/linux/$(PROBLEM_NAMES)
EXECUTABLE_NAMES := $(addprefix $(PROBLEM_NAMES), _)
EXECUTABLE_NAMES := $(addprefix $(EXECUTABLE_NAMES), ${INTERESTCALIB_DIRS})
BUILD_NUMBER_FILE := build-number.txt

interestcalib:
	$(eval i := 1)
	cd InterestCalib;\
	for dir in $(INTERESTCALIB_DIRS) ; do cd $$dir; make; cd .. ; done

interestcalib_clean:
	cd InterestCalib;\
	for dir in $(INTERESTCALIB_DIRS) ; do cd $$dir; make clean; cd .. ; done

# prints the path of the current folder.
init:
	$(info This makefile lies in: ${FINPAR_BASE_DIR})
	mkdir $(FINPAR_BASE_DIR)/bin
	mkdir $(FINPAR_BASE_DIR)/bin/linux


clean:
	rm -rf $(FINPAR_BASE_DIR)/bin/

print-%: ; @echo $* = $($*)


	# Insert something like "cp name_of_exec $(FINPAR_BASE_DIR)/bin" into above line.
	# Figure out naming convention.
	@echo $$(cat build-number.txt)
	@echo $$(($$(cat build-number.txt) + 1)) > build-number.txt
