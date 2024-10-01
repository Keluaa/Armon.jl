
CIVL ?= civl

# Add '-Xmx${CIVL_MAX_MEM}' in the Java args of the CIVL bash script for this to work
MAX_MEM ?= 64g

SOURCES = blocks.c solver.c block_interface.c time_step.c debug.c utils.c
HEADERS = blocks.h solver.h block_interface.h time_step.h debug.h utils.h vars.h
MAIN_C = main.c
MAIN_CIVL = main.cvl

MAX_ERRORS = 5

NUM_THREADS ?= 4
NUM_CYCLES ?= 3
MAX_SWEEPS ?= 2
GRID_SIZE_X ?= 3
GRID_SIZE_Y ?= 3
NPROCS ?= 4

verify_fast: export CIVL_MAX_MEM=$(MAX_MEM)
verify_fast: $(SOURCES) $(HEADERS) $(MAIN_CIVL)
	$(CIVL) verify \
		-input_mpi_nprocs=1 -input_PROC_GRID_X=1 -input_PROC_GRID_Y=1 -input_NUM_THREADS=$(NUM_THREADS) \
		-input_NUM_CYCLES=$(NUM_CYCLES) -input_MAX_SWEEPS=$(MAX_SWEEPS) \
		-input_GRID_SIZE_X=$(GRID_SIZE_X) -input_GRID_SIZE_Y=$(GRID_SIZE_Y) \
		-DSIMPLE_XCHG=1 -checkTermination \
		$(SOURCES) $(MAIN_CIVL)

replay:
	$(CIVL) replay $(SOURCES) $(MAIN_CIVL)

#-maxdepth=100000 -errorBound=$(MAX_ERRORS) \

verify_full: export CIVL_MAX_MEM=$(MAX_MEM)
verify_full: $(SOURCES) $(HEADERS) $(MAIN_CIVL)
	$(CIVL) verify \
		-input_mpi_nprocs=1 -input_PROC_GRID_X=1 -input_PROC_GRID_Y=1 \
		-DSIMPLE_XCHG=1 \
		$(SOURCES) $(MAIN_CIVL)

mpi_verify: export CIVL_MAX_MEM=$(MAX_MEM)
mpi_verify: $(SOURCES) $(HEADERS) $(MAIN_CIVL)
	$(CIVL) verify \
		-input_NUM_THREADS=$(NUM_THREADS) -input_mpi_nprocs=$(NPROCS) \
		$(SOURCES) $(MAIN_CIVL)
