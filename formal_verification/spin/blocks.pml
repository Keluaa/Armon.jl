
#include "block_grid.pml"
#include "communications.pml"

#ifndef MAX_SWEEPS
#define MAX_SWEEPS        2
#endif

#ifndef DO_TIME_STEP
#define DO_TIME_STEP      1
#endif

#ifndef DO_HALO_EXCHANGE
#define DO_HALO_EXCHANGE  1
#endif


inline init_grid()
{
    d_step {
        byte idx;
        assert(TOTAL_BLOCKS        < REMOTE_BLOCK);
        assert(TOTAL_BLOCKS        < NULL_BLOCK);
        assert(TOTAL_REMOTE_BLOCKS < REMOTE_BLOCK);
        assert(TOTAL_REMOTE_BLOCKS < NULL_BLOCK);
        assert(TOTAL_INTERFACES    < NULL_INTERFACE);

        for (idx : 0 .. (TOTAL_BLOCKS-1)) {
            block_grid[idx].state = NewCycle;
            block_grid[idx].cycle = 0;
            block_grid[idx].sweep_num = 0;
            block_grid[idx].must_wait = false;
        }

        for (idx : 0 .. (TOTAL_INTERFACES-1)) {
            block_interfaces[idx].state = XCHG_NotReady;
            block_interfaces[idx].flags = 0;
            block_interfaces[idx].is_done[0] = false;
            block_interfaces[idx].is_done[1] = false;
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

    assert(dt_contributions == TOTAL_BLOCKS);
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
        update_dt();
        new_dt_state = global_dt_state;
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


inline next_time_step(block, already_contributed)
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
        if
        :: (already_contributed) -> skip;
        :: else -> {
            // compute local_time_step
            contribute_to_dt(block);
        }
        fi
        // The first cycle requires the time step before continuing. It may have been the last block,
        // hence if DT_Done then no need to wait.
        block.must_wait = global_cycle == 0 && global_dt_state != DT_Done;
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
            :: (interface.state == XCHG_Done && interface.flags == side_flags) -> {
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


inline block_ghost_exchange(block, block_idx)
{
#if DO_HALO_EXCHANGE
    byte side;
    bool can_do_xchg, xchg_done;
    bool all_xchg_done = true;
    byte neighbour_idx[2];  // either a valid index into `block_grid` or `NULL_BLOCK`, or an index in 'remote_blocks'
    byte interface_idx[2];  // either a valid index into `block_interfaces`, `NULL_INTERFACE` or 'REMOTE_BLOCK'

    get_topology(block, block_idx);  // writes to 'neighbour_idx' and 'interface_idx'

    for (side : 0 .. 1) {
        if
#if USE_MPI
        :: (interface_idx[side] == REMOTE_BLOCK) -> {
            // 'neighbour_idx' is instead an index in 'remote_blocks'
            // TODO
        }
#endif
        :: (neighbour_idx[side] != NULL_BLOCK && interface_idx[side] != NULL_INTERFACE) -> {
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
            :: (interface_idx[side] != NULL_INTERFACE) -> { block_interfaces[interface_idx[side]].is_done[side] = false; }
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


inline block_state_machine(block, block_idx)
{
    do
    :: (block.state == NewCycle) -> {
        if
        :: (block.cycle == global_cycle) -> { block.state = TimeStep; }
        :: else -> { break; }
        fi
    }
    :: (block.state == TimeStep || block.state == InitTimeStep) -> {
        next_time_step(block, block.state == InitTimeStep);
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
        block_ghost_exchange(block, block_idx);
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
