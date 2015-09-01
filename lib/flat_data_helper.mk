# Utility makefile definitions for conveniently defining run targets
# (with some assumptions).

AD_HOC_RUNTIME_FILE=runtime.txt
AD_HOC_RESULT_FILE=result.json

AD_HOC_RUN_FLAGS=HIPERMARK_RUNTIME=$(AD_HOC_RUNTIME_FILE) HIPERMARK_RESULT=$(AD_HOC_RESULT_FILE)

run_%: $(EXECUTABLE)
	$(HIPERMARK_LIB_DIR)/linearise_data.py ../../datasets/$(subst run_,,$@)/input.json $(HIPERMARK_DATA_FIELDS) | $(AD_HOC_RUN_FLAGS) ./$(EXECUTABLE)
	@echo "Result: `cat result.json`"
	@echo "Runtime: `cat runtime.txt`ms"

run: $(EXECUTABLE)
	$(HIPERMARK_LIB_DIR)/linearise_data.py $(HIPERMARK_INPUT_FILE) $(HIPERMARK_DATA_FIELDS) | \
	 $(HIPERMARK_RUN_ENVIRONMENT) ./$(EXECUTABLE)
