
VARS ?=

MEMLIM ?= 8192
NCORES ?= 32

SWARM_PROC ?= 128
SWARM_TIME ?= 1.0
SWARM_MEM ?= 1024M

# TODO: remove channel assertions validity checks just to be sure for now, re-enable if the performance gain is noticeable
# XUSAFE ?= -DXUSAFE
XUSAFE ?=

SOURCES = utils.pml blocks.pml block_grid.pml solver.pml main.pml
GRID_DEF = ./grid_definition.c ./grid_definition.h ./grid_definition.pml

$(GRID_DEF): ./gen_grid.jl
	@echo "use 'julia ./gen_grid.jl' to (re)generate the topology of the block grid" && false

./pan.c: $(SOURCES) $(GRID_DEF)
	spin -a $(VARS) ./main.pml

./pan_exhaustive: ./pan.c
	$(CC) -O2 -DMEMLIM=65536 -DVECTORSZ=4160 $(XUSAFE) -DSAFETY -DNOCLAIM -DCOLLAPSE -DREACH -w -o pan_exhaustive ./pan.c

./pan_bitstate: ./pan.c
	$(CC) -O2 -DMEMLIM=$(MEMLIM) -DBITSTATE -DVECTORSZ=4160 $(XUSAFE) -DNOCLAIM -w -o pan_bitstate ./pan.c

./pan_bitstate_non_progress: ./pan.c
	$(CC) -O2 -DMEMLIM=$(MEMLIM) -DBITSTATE -DVECTORSZ=4160 $(XUSAFE) -DNP -w -o pan_bitstate_non_progress ./pan.c

./multi_core:
	mkdir -p ./multi_core

./multi_core/pan.c: $(SOURCES) $(GRID_DEF) ./multi_core
	 cd ./multi_core && spin -a -DNO_HIDDEN $(VARS) ../main.pml

./multi_core/pan_exhaustive: ./multi_core/pan.c ./multi_core
	cd ./multi_core && $(CC) -O2 -DMEMLIM=$(MEMLIM) -DVECTORSZ=4160 $(XUSAFE) -DSAFETY -DNOCLAIM -DCOLLAPSE -DREACH -DNCORE=$(NCORES) -DPMAX=5 -DVMAX=344 -w -o pan_exhaustive ./pan.c

./multi_core/pan_bitstate: ./multi_core/pan.c ./multi_core
	cd ./multi_core && $(CC) -O2 -DMEMLIM=$(MEMLIM) -DBITSTATE -DVECTORSZ=4160 $(XUSAFE) -DNOCLAIM -DNCORE=$(NCORES) -DPMAX=5 -DVMAX=344 -w -o pan_bitstate ./pan.c

./swarm_check:
	mkdir -p ./swarm_check
	mkdir -p ./swarm_check/trails

./swarm_check/swarm_check.swarm: export CCOMMON=-O2 -DVECTORSZ=4160 $(XUSAFE) -DNOCLAIM
./swarm_check/swarm_check.swarm: export RCOMMON=
./swarm_check/swarm_check.swarm: ./swarm_check
	# Bu default, 128 processes running for 1h on 128 cores using up to 8GB each and exploring about 500000 states per second
	swarm -c$(SWARM_PROC) -m$(SWARM_MEM) -t$(SWARM_TIME) -s500000 -d10000000 -f swarm_check
	mv ./swarm_check.swarm ./swarm_check/swarm_check.swarm

single_grid:
	julia ./gen_grid.jl ./grids/default_grid_3x3_4threads.toml

mpi_grid:
	julia ./gen_grid.jl ./grids/proc=2x2_grid=3x3_4threads.toml

preprocessor: $(GRID_DEF)
	spin -I ./main.pml

pan: ./pan.c

pan_multi_core: ./multi_core/pan.c

check: ./pan_exhaustive
	./pan_exhaustive -m10000000 -c1 -q -I

check_approx: ./pan_bitstate
	./pan_bitstate -m10000000 -c1 -k1 -w34 -q -I

check_non_progress: ./pan_bitstate_non_progress
	./pan_bitstate_non_progress -m10000000 -l -c1 -k1 -w34 -q -I

check_multi_core: ./multi_core/pan_exhaustive
	cd ./multi_core && ./pan_exhaustive -m10000000 -c1 -q -I

approx_check_multi_core: ./multi_core/pan_bitstate
	cd ./multi_core && ./pan_bitstate -m10000000 -c1 -w30 -q

full_trail: ./main.pml.trail
	./pan_exhaustive -q -r ./main.pml.trail

approx_trail: ./main.pml.trail
	./pan_bitstate -q -r ./main.pml.trail

non_progress_trail: ./main.pml.trail
	./pan_bitstate_non_progress -l -q -r ./main.pml.trail

swarm: ./pan.c ./swarm_check/swarm_check.swarm
	cp -f ./pan.* ./swarm_check/
	rm -f ./swarm_check/*.trail
	cd ./swarm_check && ./swarm_check.swarm
	if [ -f ./swarm_check/*.trail ]; then mv ./swarm_check/*.trail ./swarm_check/trails/; fi
	@echo "Found $(shell ls -1 ./swarm_check/trails | wc -l) runs with problems. Trails can be found in './swarm_check/trails'"
	$(shell ls -1 ./swarm_check/trails | [[ $$(wc -l) -eq 0 ]])

swarm_trail: ./pan.c ./swarm_check/swarm_check.swarm
	cp -f ./pan.* ./swarm_check/
	cd ./swarm_check && ./swarm_check.swarm prep
	@echo "Use './swarm_check/pan[1-4] -r ./swarm_check/trails/main.pml_pan[1-4]_<seed>.trail' to follow a trail"

clean_swarm:
	rm -rf ./swarm_check

clean_multi_core:
	rm -rf ./multi_core

clean_pan:
	rm -f ./pan* ./*.trail

clean: clean_pan clean_swarm clean_multi_core
	rm -f $(GRID_DEF)

help:
	@echo "Formal verification of the solver for multithreading + MPI using SPIN"
	@echo "SPIN can be downloaded by following those instructions: https://spinroot.com/spin/Man/README.html"
	@echo "Targets:"
	@echo "   single_grid            Create './grid_definition.c' for a single grid of 3x3 blocks using 4 threads"
	@echo "   mpi_grid               Create './grid_definition.c' for 4 MPI ranks each with 3x3 blocks and 4 threads"
	@echo "   preprocessor           Output the Pomela source after the preprocessor pass"
	@echo "   pan                    Compile the Pomela source to C"
	@echo "   check                  Perform an exhaustive simulation of all solver states (impossibly slow)"
	@echo "   check_approx           Perform an approximate simulation of some solver states (slow)"
	@echo "   check_non_progress     Perform an approximate simulation, searching for non-progress cycles (slow)"
	@echo "   swarm                  Perform an approximate simulation by using multiple processes (fast)"
	@echo "   trail                  Print the trail and final state of an exhaustive simulation which failed"
	@echo "   approx_trail           Print the trail and final state of an approximate simulation which failed"
	@echo "   non_progress_trail     Print the trail and final state of an approximate simulation for non-progress cycles which failed"
	@echo "   swarm_trail            Setup necessary files to follow the trails generated during a swarm run"
	@echo "   help                   Display this message"
	@echo "   clean                  Remove all generated files"
	@echo "Verification variables (specify with 'make VARS=\"-DVAR_NAME=value\" pan'):"
	@echo "   NUM_CYCLES        [3]  Number of cycles to do"
	@echo "   MAX_SWEEPS        [2]  Number of block sweeps to do per cycle"
	@echo "   DO_TIME_STEP      [1]  If the reduction for time step at each cycle should be simulated"
	@echo "   DO_HALO_EXCHANGE  [1]  If the halo exchange of each block for each sweep should be simulated"
	@echo "   USE_MPI           [0]  If MPI communications should be simulated"
	@echo "Verification ('check*' targets) compilation options:"
	@echo "   MEMLIM         [8192]  Maximum amount of memory to use, in Megabytes"
	@echo "Swarm options (specify with 'make SWARM_*=value swarm'):"
	@echo "   SWARM_PROC      [128]  Number of processes to use"
	@echo "   SWARM_TIME      [1.0]  Maximum time for the verification, in hours"
	@echo "   SWARM_MEM     [1024M]  Maximum amount of memory available for each process"
