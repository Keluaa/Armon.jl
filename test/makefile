
PROJECT_DIR = ..
TEST_DIR = .
DATA_REF_DIR = $(TEST_DIR)/reference_data

MPI_PROCS ?= 25
MPI_CORES ?= 4
MPI_EXEC ?= mpiexec
MPI_EXEC_ARGS ?= -n $(MPI_PROCS) -c $(MPI_CORES)
MPI_JL_ARGS ?= --color=yes

JULIA_THREADS ?= $(MPI_CORES)
JULIA ?= julia
JULIA_ARGS ?= -O3 --project=$(PROJECT_DIR) -t $(JULIA_THREADS)

main_test_file = $(TEST_DIR)/runtests.jl
create_data_ref_file = $(DATA_REF_DIR)/create_references.jl


# When target is 'test', forward other args to Julia
ifeq (test,$(firstword $(MAKECMDGOALS)))
  TEST_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(TEST_ARGS):;@:)
endif


.PHONY: create_solver_reference short perf mpi all

create_solver_reference:
	$(JULIA) $(JULIA_ARGS) $(create_data_ref_file)

mpi:
	$(MPI_EXEC) $(MPI_EXEC_ARGS) $(JULIA) $(JULIA_ARGS) $(MPI_JL_ARGS) $(main_test_file) mpi

test:
	$(JULIA) $(JULIA_ARGS) $(main_test_file) $(TEST_ARGS)

all:
	$(JULIA) $(JULIA_ARGS) $(main_test_file) all
