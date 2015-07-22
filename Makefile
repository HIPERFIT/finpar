FINPAR_BASE_DIR := $(shell pwd)
INTERESTCALIB_BIN_DIR := $(FINPAR_BASE_DIR)/bin/linux/InterestCalib
INTERESTCALIB_DIRS := CppOpenMP/ OrigCpp/

interestcalib: 
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
	mkdir $(INTERESTCALIB_BIN_DIR)


clean:
	rm -rf $(FINPAR_BASE_DIR)/bin/
	cd Interest

print-%: ; @echo $* = $($*)


	# Insert something like "cp name_of_exec $(FINPAR_BASE_DIR)/bin" into above line.
	# Figure out naming convention.
