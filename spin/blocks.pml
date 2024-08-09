
#ifndef MAX_SWEEPS
#define MAX_SWEEPS    2
#endif

#ifndef GRID_SIZE_X
#define GRID_SIZE_X   3
#endif
#ifndef GRID_SIZE_Y
#define GRID_SIZE_Y   3
#endif

#define TOTAL_BLOCKS         GRID_SIZE_X*GRID_SIZE_Y
#define TOTAL_INTERFACES     2*TOTAL_BLOCKS-GRID_SIZE_X-GRID_SIZE_Y
#define INTERFACES_Y_OFFSET  TOTAL_BLOCKS-GRID_SIZE_X-1
#define NO_NEIGHBOUR         255
#define NO_INTERFACE         32767

#ifndef DO_TIME_STEP
#define DO_TIME_STEP      1
#endif
#ifndef DO_HALO_EXCHANGE
#define DO_HALO_EXCHANGE  1
#endif

#ifndef USE_MPI
#define USE_MPI           0
#endif
#ifndef MPI_CHANNEL_SIZE
#define MPI_CHANNEL_SIZE  10
#endif


mtype:BlockStates = { NewCycle, TimeStep, InitTimeStep, NewSweep, EOS, Exchange, Fluxes, CellUpdate, Remap, EndCycle };
mtype:BlockXChg   = { XCHG_NotReady, XCHG_InProgress, XCHG_Done };

typedef Block {
    mtype:BlockStates state;
    byte cycle;  // local cycle of the block
    byte sweep_num;
    bool must_wait;  // if the block is waiting for other blocks
    byte pos[2];
};

typedef BlockInterface {
    // in the Julia implementation, 'state' and 'flags' are stored in the same byte,
    // therefore atomic operations on both values at the same time is possible.
    mtype:BlockXChg state;
    byte flags;
    bool is_done[2];
};

Block block_grid[TOTAL_BLOCKS];
BlockInterface block_interfaces[TOTAL_INTERFACES]

// MPI neighbours
chan subdomain_neighbours_dt_reduction     = [MPI_CHANNEL_SIZE] of { byte };
chan subdomain_neighbours_halo_exchange[4] = [MPI_CHANNEL_SIZE] of { byte };

// Global state
mtype:TimeStepState = { DT_Ready, DT_AllContributed, DT_DoingMPI, DT_WaitingForMPI, DT_Done };
mtype:TimeStepState global_dt_state = DT_Ready;
byte dt_contributions = 0;  // number of blocks which contributed to the time step calculation for this cycle
byte global_cycle = 0;  // the current cycle of the whole solver


inline print_block(block)
{
    printf("block at (%d,%d), state=", block.pos[0], block.pos[1]);
    printm(block.state);
    printf(", cycle=%d, sweep=%d", block.cycle, block.sweep_num);
}


inline init_grid()
{
    byte i, j, idx;
    assert(TOTAL_BLOCKS < NO_NEIGHBOUR);
    assert(TOTAL_INTERFACES < NO_INTERFACE);
    d_step {
        for (j : 0 .. (GRID_SIZE_Y-1)) {
            for (i : 0 .. (GRID_SIZE_X-1)) {
                idx = j * GRID_SIZE_X + i;
                block_grid[idx].state = NewCycle;
                block_grid[idx].cycle = 0;
                block_grid[idx].sweep_num = 0;
                block_grid[idx].must_wait = false;
                block_grid[idx].pos[0] = i;
                block_grid[idx].pos[1] = j;
            }
        }

        for (i : 0 .. (TOTAL_INTERFACES-1)) {
            block_interfaces[i].state = XCHG_NotReady;
            block_interfaces[i].flags = 0;
            block_interfaces[i].is_done[0] = false;
            block_interfaces[i].is_done[1] = false;
        }
    };
}


inline update_dt()
{
    mtype:TimeStepState prev_dt_state = global_dt_state;
    if
    :: (prev_dt_state == DT_AllContributed) -> {
#if USE_MPI
        // TODO: start the `MPI_Iallreduce`
        global_dt_state = DT_DoingMPI;
        goto skip_dt_update;
#else
        skip;
#endif
    }
    :: (prev_dt_state == DT_WaitingForMPI) -> {
        // receive the MPI reduction result
        assert(nempty(subdomain_neighbours_dt_reduction));
        subdomain_neighbours_dt_reduction?_;
    }
    :: else -> assert(false);
    fi

    dt_contributions = 0;
    global_dt_state = DT_Done;

skip_dt_update:
    skip;
}


inline wait_for_dt(new_dt_state)
{
    bool cas_ok = false;
    atomic_cas(cas_ok, global_dt_state, DT_DoingMPI, DT_WaitingForMPI);
    if
    :: (cas_ok) -> {
#if USE_MPI
        // TODO: wait until `subdomain_neighbours_dt_reduction` has a value + make sure only one thread waits on the request
#endif
        update_dt()
        new_dt_state = global_dt_state;  // TODO: is this correct? shouldn't we use a "return value" from `update_dt` instead?
    }
    :: else -> { new_dt_state = DT_WaitingForMPI; };
    fi
}


inline contribute_to_dt(block)
{
    byte current_contributions = 0;
    bool cas_ok = false;
    atomic {
        dt_contributions++;
        current_contributions = dt_contributions;
    };

    if
    :: (current_contributions == TOTAL_BLOCKS) -> {
        atomic_cas(cas_ok, global_dt_state, DT_Ready, DT_AllContributed);
        if
        :: (cas_ok) -> update_dt();
        :: else -> skip;
        fi
    }
    :: else -> skip;
    fi
}


inline next_cycle()
{
    mtype:TimeStepState current_dt_state;

#if DO_TIME_STEP
    current_dt_state = global_dt_state;
retry_next_dt:
    if
    :: (current_dt_state == DT_DoingMPI) -> {
        wait_for_dt(current_dt_state);
        goto retry_next_dt;
    }
    :: (current_dt_state == DT_Done) -> {
        global_dt_state = DT_Ready;
    }
    :: else -> assert(false);  // the global time step must be done before continuing
    fi
#endif

    global_cycle++;
}


inline next_time_step(block)
{
#if DO_TIME_STEP
    mtype:TimeStepState dt_state = global_dt_state;
retry_time_step:
    if
    :: (dt_state == DT_DoingMPI) -> {
        wait_for_dt(dt_state);
        assert(dt_state != DT_DoingMPI);
        goto retry_time_step;
    }
    :: (dt_state == DT_Ready) -> {
        // local_time_step
        contribute_to_dt(block);
        block.must_wait = global_cycle == 0;  // The first cycle requires the time step before continuing
    }
    :: (dt_state == DT_Done) -> {
        block.must_wait = false;
    }
    :: else -> { block.must_wait = true; }
    fi
#else
    block.must_wait = false;
#endif
}


inline mark_ready_for_exchange(interface, side, can_do_xchg, xchg_done)
{
    mtype:BlockXChg int_state;
    byte int_flags;
    byte side_flags = (side == 0 -> 2 : 1);
    atomic {
        int_state = interface.state;
        int_flags = interface.flags;
    };

    if
    :: (int_state == XCHG_InProgress) -> {
        can_do_xchg = false;
        xchg_done = false;
        goto xchg_marked;
    }
    :: (int_state == XCHG_Done) -> {
        can_do_xchg = false;
        // interface_acknowledge_exchange (a CAS but also reset the flags on success)
        atomic {
            if
            :: (interface.state == XCHG_Done && interface.flags != side_flags) -> {
                interface.state = XCHG_NotReady;
                interface.flags = 0;
                xchg_done = true;
            }
            :: else -> { xchg_done = false; }
            fi
        }
        goto xchg_marked;
    }
    :: (int_state == XCHG_NotReady) -> {
        can_do_xchg = false;
        xchg_done = false;
    };
    :: else -> assert(false);
    fi

    if
    :: ((int_flags & side_flags) == 0) -> {
        // interface_side_ready
        atomic {
            interface.flags = interface.flags | int_flags | side_flags;
            int_state = interface.state;
            int_flags = interface.flags;
        }
    };
    :: else -> skip;
    fi

    if
    :: (int_flags == 3) -> {
        // interface_start_exchange (CAS with the state and flags of the interface)
        atomic {
            if
            :: (interface.state == XCHG_NotReady && interface.flags == 3) -> {
                interface.state = XCHG_InProgress;
                interface.flags = (side == 0 -> 1 : 2);  // opposite of side_flags
                can_do_xchg = true;
            }
            :: else -> { can_do_xchg = false; }
            fi
        }
    }
    :: else -> { can_do_xchg = false; };
    fi

xchg_marked:
    skip;
}


inline block_ghost_exchange(block)
{
#if DO_HALO_EXCHANGE
    byte side;
    bool can_do_xchg, xchg_done;
    bool all_xchg_done = true;
    short interface_idx[2];  // either a valid index into `block_interfaces` or `NO_INTERFACE`
    byte  neighbour_idx[2];  // either a valid index into `block_grid` or `NO_NEIGHBOUR`
    d_step {
        if
        :: (block.sweep_num == 0) -> {
            // Halo exchange along the X axis
            neighbour_idx[0] = ( (block.pos[0] > 0)             -> (block.pos[1] * GRID_SIZE_X + block.pos[0] - 1) : NO_NEIGHBOUR );
            neighbour_idx[1] = ( (block.pos[0] < GRID_SIZE_X-1) -> (block.pos[1] * GRID_SIZE_X + block.pos[0] + 1) : NO_NEIGHBOUR );
            interface_idx[0] = ( (block.pos[0] > 0)             -> (block.pos[1] * (GRID_SIZE_X-1) + block.pos[0] - 1) : NO_INTERFACE );
            interface_idx[1] = ( (block.pos[0] < GRID_SIZE_X-1) -> (block.pos[1] * (GRID_SIZE_X-1) + block.pos[0] + 1) : NO_INTERFACE );
        }
        :: (block.sweep_num == 1) -> {
            // Halo exchange along the Y axis
            neighbour_idx[0] = ( (block.pos[1] > 0)             -> ((block.pos[1]-1) * GRID_SIZE_X + block.pos[0]) : NO_NEIGHBOUR );
            neighbour_idx[1] = ( (block.pos[1] < GRID_SIZE_Y-1) -> ((block.pos[1]+1) * GRID_SIZE_X + block.pos[0]) : NO_NEIGHBOUR );
            interface_idx[0] = ( (block.pos[1] > 0)             -> (INTERFACES_Y_OFFSET + block.pos[0] * (GRID_SIZE_Y-1) + block.pos[1] - 1) : NO_INTERFACE );
            interface_idx[1] = ( (block.pos[1] < GRID_SIZE_Y-1) -> (INTERFACES_Y_OFFSET + block.pos[0] * (GRID_SIZE_Y-1) + block.pos[1] + 1) : NO_INTERFACE );
        }
        :: else -> assert(false);
        fi

        // Having no neighbour on one side must imply that there is no interface on that side
        printf("block at (%d, %d) and sweep %d will xchg with:\n - left:  idx=%d, int=%d\n - right: idx=%d, int=%d\nthis is no int: %d and offset: %d\n",
               block.pos[0], block.pos[1], block.sweep_num, neighbour_idx[0], interface_idx[0], neighbour_idx[1], interface_idx[1], NO_INTERFACE, INTERFACES_Y_OFFSET);
        assert(((neighbour_idx[0] == NO_NEIGHBOUR) ^ (interface_idx[0] == NO_INTERFACE)) == 0);
        assert(((neighbour_idx[1] == NO_NEIGHBOUR) ^ (interface_idx[1] == NO_INTERFACE)) == 0);
    };

    for (side : 0 .. 1) {
        if
        :: (neighbour_idx[side] != NO_NEIGHBOUR && interface_idx[side] != NO_INTERFACE) -> {
            if
            :: (block_interfaces[interface_idx[side]].is_done[side]) -> skip;  // side is already done
            :: else -> {
                // block_ghost_exchange between two local blocks
                mark_ready_for_exchange(block_interfaces[interface_idx[side]], side, can_do_xchg, xchg_done);
                if
                ::(can_do_xchg) -> {
                    // "do the exchange between the blocks"
                    // exchange_done
                    atomic {
                        block_interfaces[interface_idx[side]].state = XCHG_Done;
                    }
                    xchg_done = true;
                }
                :: else -> skip;
                fi

                if
                :: (xchg_done) -> { block_interfaces[interface_idx[side]].is_done[side] = true; }
                :: else -> { all_xchg_done = false; }
                fi
            }
            fi
        }
        :: else -> skip;
        fi
    }

    if
    :: (all_xchg_done) -> {
        // Reset the interfaces
        for (side : 0 .. 1) {
            if
            :: (interface_idx[side] != NO_INTERFACE) -> { block_interfaces[interface_idx[side]].is_done[side] = false; }
            :: else -> skip;
            fi
        }
    }
    :: else -> skip;
    fi

    block.must_wait = !all_xchg_done;
#else
    block.must_wait = false;
#endif
}


inline block_state_machine(block)
{
    do
    :: (block.state == NewCycle) -> {
        if
        :: (block.cycle == global_cycle) -> { block.state = TimeStep; }
        :: else -> { break; }
        fi
    }
    :: (block.state == TimeStep || block.state == InitTimeStep) -> {
        next_time_step(block);
        if
        :: (block.must_wait) -> { block.state = InitTimeStep; break; }
        :: else -> { block.state = NewSweep; }
        fi
    }
    :: (block.state == NewSweep) -> {
        if
        :: (block.sweep_num == MAX_SWEEPS-1) -> { block.sweep_num = 0; block.state = EndCycle; }
        :: else -> { block.sweep_num++; block.state = EOS; }
        fi
    }
    :: (block.state == EOS) -> {
        block.state = Exchange;
    }
    :: (block.state == Exchange) -> {
        block_ghost_exchange(block);
        if
        :: (block.must_wait) -> break;
        :: else -> { block.state = Fluxes; }
        fi
    }
    :: (block.state == Fluxes)     -> { block.state = CellUpdate; }
    :: (block.state == CellUpdate) -> { block.state = Remap;      }
    :: (block.state == Remap)      -> { block.state = NewSweep;   }
    :: (block.state == EndCycle) -> {
        block.cycle++;
        block.state = NewCycle;
        break;
    }
    :: else -> assert(false);
    od
}
