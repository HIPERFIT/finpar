# This makefile provides utility definitions for having either
# REAL_IS_FLOAT or REAL_IS_DOUBLE defined during compilation.  This is
# based on the environment variable HIPERMARK_CONFIG_REAL_TYPE, which
# must be defined as a static configuration option.  If
# HIPERMARK_CONFIG_REAL_TYPE is not set at all, it will be assumed to
# be 'double'
#
# This makefile fragment appends to CXXFLAGS.

ifeq ($(HIPERMARK_CONFIG_REAL_TYPE),float)
	CXXFLAGS += -DREAL_IS_FLOAT
else ifeq ($(HIPERMARK_CONFIG_REAL_TYPE),double)
	CXXFLAGS += -DREAL_IS_DOUBLE
else
	CXXFLAGS += -DREAL_IS_DOUBLE
endif
