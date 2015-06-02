# can be set using $>gmake GHC=/path/to/my/ghc 
ifeq ($(strip $(GHC)),) 
  GHC=`which ghc`
endif

# switch this on if llvm is enabled
#BACKEND=-fllvm

OPTIMIZATION=-O2 -msse2 -rtsopts 
#-Wall
EXECUTABLE=VolCalib

all: $(EXECUTABLE)

$(EXECUTABLE): % : %.hs
	$(GHC) $(OPTIMIZATION) --make $< -o $@

run_small: $(EXECUTABLE)
	cat ../Data/Small/input.data ../Data/Small/output.data | ./$(EXECUTABLE)

run_medium:  $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE)

# the Large dataset runs out of memory on the current Haskell implementation!
run_large:  $(EXECUTABLE)
	cat ../Data/Medium/input.data ../Data/Medium/output.data | ./$(EXECUTABLE)
#	cat ../Data/Large/input.data ../Data/Large/output.data | ./$(EXECUTABLE)


#-------------------------------------------------------------------------
# clean
#-------------------------------------------------------------------------
.PHONY  : clean_executables clean_obj clean_source_backup
.SILENT : clean_executables clean_obj clean_source_backup

clean: clean_executables clean_obj clean_source_backup

clean_executables:
	rm -f $(EXECUTABLE)

clean_obj :
	rm -f *.o *.hi

clean_source_backup:
	rm -f *~

