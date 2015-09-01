# Permit standalone compilation.
HIPERMARK_LIB_DIR ?= ../../../../lib/

# Default to double if nothing else defined.
HIPERMARK_CONFIG_REAL_TYPE ?= double

include $(HIPERMARK_LIB_DIR)/setup.mk
include $(HIPERMARK_LIB_DIR)/setup_real_type.mk

# FIXME - this should not even exist!
include $(HIPERMARK_LIB_DIR)/platform.mk

INCLUDES    += -I$(HIPERMARK_BENCHMARK_LIB_DIR)/include
INCLUDES    += -I$(HIPERMARK_LIB_DIR)/include

GPU_OPTS   = -D lgWARP=$(GPU_LG_WARP) -D GPU_DEV_ID=$(GPU_DEVICE_ID) \
	     -D GPU_LOC_MEM=$(GPU_LOCAL_MEM) -D GPU_CONST_MEM=$(GPU_CONST_MEM) \
	     -D GPU_REG_MEM=$(GPU_REG_MEM) -D GPU_GLB_MEM=$(GPU_DEVICE_MEM) \
	     -D GPU_TILE=$(GPU_LOCAL_MEM_PER_TH) -D GPU_CORES=$(GPU_NUM_CORES) \
			 -D CURR_DIR_PATH='"$(MAKE_DIR)"'

SOURCES_CPP =$(HIPERMARK_IMPLEMENTATION_DIR)/VolCalibAll.cpp
HELPERS     =PrepareKernels.h Vect_CPU.h Vect_GPU.h VolCalibInit.h \
	     $(HIPERMARK_BENCHMARK_LIB_DIR)/Constants.h \
	     $(HIPERMARK_BENCHMARK_LIB_DIR)/ParseInput.h \
	     $(HIPERMARK_BENCHMARK_LIB_DIR)/DataStructConst.h \
	     $(HIPERMARK_BENCHMARK_LIB_DIR)/ParPrefixUtil.h
EXECUTABLE  =VolCalibAll

$(EXECUTABLE):
	cp $(HIPERMARK_IMPLEMENTATION_DIR)/*cl .
	$(CXX) $(CXXFLAGS) -D IS_GPU=1 $(GPU_OPTS) $(INCLUDES) -o $(EXECUTABLE) $(SOURCES_CPP) $(LIB)


clean:
	rm -f $(EXECUTABLE) $(OBJECTS) ./*.ptx ./Debug.txt
	@# clean nVidia compiler cache
	rm -rf $(HOME)/.nv/ComputeCache/*

HIPERMARK_DATA_FIELDS=outer num_x num_y num_t
ifdef HIPERMARK_CONFIG_NUM_THREADS
  HIPERMARK_RUN_ENVIRONMENT=OMP_NUM_THREADS=$(HIPERMARK_CONFIG_NUM_THREADS)
endif
include $(HIPERMARK_LIB_DIR)/flat_data_helper.mk
