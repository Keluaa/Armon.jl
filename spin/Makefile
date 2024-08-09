
VARS ?=

SOURCES = utils.pml blocks.pml threads_workload.pml solver.pml

./pan.c: $(SOURCES)
	spin -a $(VARS) ./solver.pml

./pan_exhaustive: ./pan.c
	$(CC) -DMEMLIM=4096 -O2 -DXUSAFE -DSAFETY -DNOCLAIM -DCOLLAPSE -w -o pan_exhaustive ./pan.c

./pan_bitstate: ./pan.c
	$(CC) -DMEMLIM=4096 -O2 -DBITSTATE -DXUSAFE -DNP -DNOCLAIM -w -o pan_bitstate ./pan.c

preprocessor:
	spin -I ./solver.pml

pan: ./pan.c

check: ./pan_exhaustive
	./pan_exhaustive

approx_check: ./pan_bitstate
	./pan_bitstate -m10000000 -k5 -l -c1 -w34

follow_full_trail: ./pan_exhaustive ./solver.pml.trail
	./pan_exhaustive -r ./solver.pml.trail

follow_approx_trail: ./pan_bitstate ./solver.pml.trail
	./pan_bitstate -r ./solver.pml.trail

clean:
	rm -f ./pan* ./*.trail

help:
	@echo "Formal verification of the solver for multithreading + MPI using SPIN"
	@echo "SPIN can be downloaded by following those instructions: https://spinroot.com/spin/Man/README.html"
	@echo "Targets:"
	@echo "   preprocessor           Output the Pomela source after the preprocessor pass"
	@echo "   pan                    Compile the Pomela source to C"
	@echo "   check                  Perform an exhaustive simulation of all solver states (slow)"
	@echo "   approx_check           Perform an approximate simulation of most solver states (fast)"
	@echo "   follow_trail           Print the trail and final state of an exhaustive simulation which failed"
	@echo "   follow_approx_trail    Print the trail and final state of an approximate simulation which failed"
	@echo "   help                   Display this message"
	@echo "   clean                  Remove all generated files"
	@echo "Verification variables (specify with \"VARS='-DVAR_NAME=value'\"):"
	@echo "   NUM_THREADS       [4]  Number of threads to simulate"
	@echo "   NUM_CYCLES        [3]  Number of cycles to do"
	@echo "   MAX_SWEEPS        [2]  Number of block sweeps to do per cycle"
	@echo "   GRID_SIZE_X       [3]  Size of the block grid along X"
	@echo "   GRID_SIZE_Y       [3]  Size of the block grid along Y"
	@echo "                          GRID_SIZE_X*GRID_SIZE_Y should not be >= 255"
	@echo "   DO_TIME_STEP      [1]  If the reduction for time step at each cycle should be simulated"
	@echo "   DO_HALO_EXCHANGE  [1]  If the halo exchange of each block for each sweep should be simulated"
	@echo "   USE_MPI           [0]  If MPI communications should be simulated"
	@echo "   CHECK_DISTRIB     [0]  Perform checks ensuring all blocks are correctly distributed among all threads"
