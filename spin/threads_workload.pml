
#ifndef CHECK_DISTRIB
#define CHECK_DISTRIB 0
#endif


typedef ThreadWorkload {
    byte num_blocks;  // number of blocks assigned to the thread
    byte blocks[TOTAL_BLOCKS];  // `0` to `TOTAL_BLOCKS-1` for blocks assigned to the thread, `NO_NEIGHBOUR` otherwise
};
ThreadWorkload threads_workload[NUM_THREADS];


#if CHECK_DISTRIB
hidden byte all_blocks[TOTAL_BLOCKS];
byte bad = 0;
#endif
inline init_threads_workload()
{
    d_step {
        byte blocks_per_thread = TOTAL_BLOCKS / NUM_THREADS;
        byte remaining_blocks = TOTAL_BLOCKS - NUM_THREADS * blocks_per_thread;
        byte prev_tids_blocks, tid_blocks;
        byte _idx, _tid, _i_idx;

        // Distribute all blocks to all threads evenly
        for (_tid : 0 .. NUM_THREADS-1) {
            prev_tids_blocks = blocks_per_thread * _tid;
            tid_blocks = blocks_per_thread;
            if
            :: (_tid >= remaining_blocks) -> { prev_tids_blocks = prev_tids_blocks + remaining_blocks; }
            :: else -> { prev_tids_blocks = prev_tids_blocks + _tid; tid_blocks++; }
            fi

            threads_workload[_tid].num_blocks = tid_blocks;
            for (_idx : 0 .. tid_blocks-1) {
                threads_workload[_tid].blocks[_idx] = _idx + prev_tids_blocks;
            }
            for (_idx : tid_blocks .. TOTAL_BLOCKS-1) {
                threads_workload[_tid].blocks[_idx] = NO_NEIGHBOUR;
            }
        }

#if CHECK_DISTRIB
        // Check that all blocks are assigned to a single thread

        for (_idx : 0 .. TOTAL_BLOCKS-1) {
            all_blocks[_idx] = NO_NEIGHBOUR;
        }

        for (_tid : 0 .. NUM_THREADS-1) {
            for (_i_idx : 0 .. threads_workload[_tid].num_blocks-1) {
                _idx = threads_workload[_tid].blocks[_i_idx];
                if
                :: (_idx == NO_NEIGHBOUR) -> skip;
                :: (all_blocks[_idx] == NO_NEIGHBOUR) -> { all_blocks[_idx] = _tid; }
                :: else -> {
                    printf("Workload: _idx %d is shared by threads %d and %d\n", _idx, all_blocks[_idx], _tid);
                    bad++;
                }
                fi
            }
        }

        for (_idx : 0 .. TOTAL_BLOCKS-1) {
            if
            :: (all_blocks[_idx] == NO_NEIGHBOUR) -> {
                printf("Workload: _idx %d is not assigned to a thread\n", _idx, all_blocks[_idx]);
                bad++;
            }
            :: else -> skip;
            fi
        }
#endif
    }

#if CHECK_DISTRIB
    assert(bad == 0);
#endif
}
