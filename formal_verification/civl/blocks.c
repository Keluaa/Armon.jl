
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "blocks.h"
#include "time_step.h"
#include "debug.h"
#include "mpi_comms.h"

AtomicVar enum TimeStepState global_dt_state = DT_Ready;
AtomicVar byte dt_contributions = 0;
byte global_cycle = 0;

static void init_thread_workload(struct ThreadWorkload* thread_workload, struct Block* all_blocks, int num_threads, int tid)
{
    int blocks_per_thread = TOTAL_BLOCKS / num_threads;
    int remaining_blocks = TOTAL_BLOCKS - num_threads * blocks_per_thread;

    int prev_tids_blocks = blocks_per_thread * tid;
    int tid_blocks = blocks_per_thread;
    if (tid >= remaining_blocks) {
        prev_tids_blocks += remaining_blocks;
    } else {
        prev_tids_blocks += tid;
        tid_blocks++;
    }

    thread_workload->num_blocks = tid_blocks;
    thread_workload->threads_blocks = malloc(sizeof(struct Block*) * tid_blocks);
    for (int block_i = 0; block_i < tid_blocks; block_i++) {
        thread_workload->threads_blocks[block_i] = &all_blocks[prev_tids_blocks + block_i];
    }
}

struct BlockGrid* init_grid(MPI_Comm comm, const int* rank_pos, const int* neighbour_ranks, int num_threads)
{
    struct BlockGrid* block_grid = malloc(sizeof(struct BlockGrid));
    block_grid->blocks = malloc(sizeof(struct Block) * TOTAL_BLOCKS);
    block_grid->remote_blocks = malloc(sizeof(struct RemoteBlock) * TOTAL_REMOTE_BLOCKS);
    block_grid->interfaces = malloc(sizeof(struct BlockInterface) * TOTAL_INTERFACES);

    block_grid->global_origin[0] = rank_pos[0] * GRID_SIZE_X;
    block_grid->global_origin[1] = rank_pos[1] * GRID_SIZE_Y;

    for (int i = 0; i < TOTAL_INTERFACES; i++) {
        struct BlockInterface* interface = &block_grid->interfaces[i];
#if SIMPLE_XCHG
        interface->ready[0] = false;
        interface->ready[1] = false;
        interface->bint_state = XCHG_NotReady;
#else
        interface->int_state = 0;
#endif
        interface->is_done[0] = false;
        interface->is_done[1] = false;
    }

    const int interfaces_X_offset = 0;
    const int interfaces_Y_offset = TOTAL_BLOCKS - GRID_SIZE_X;
    int global_grid_size[] = { PROC_GRID_X * GRID_SIZE_X, PROC_GRID_Y * GRID_SIZE_Y };
    int remote_block_i = 0;
    for (byte j = 0; j < GRID_SIZE_Y; j++) {
        for (byte i = 0; i < GRID_SIZE_X; i++) {
            byte idx = j * GRID_SIZE_X + i;

            struct BlockInterface *left_int, *right_int, *bottom_int, *top_int;
            left_int   = (i > 0)             ? &block_grid->interfaces[interfaces_X_offset + j * (GRID_SIZE_X-1) + i - 1] : NULL;
            right_int  = (i < GRID_SIZE_X-1) ? &block_grid->interfaces[interfaces_X_offset + j * (GRID_SIZE_X-1) + i    ] : NULL;
            bottom_int = (j > 0)             ? &block_grid->interfaces[interfaces_Y_offset + i * (GRID_SIZE_Y-1) + j - 1] : NULL;
            top_int    = (j < GRID_SIZE_Y-1) ? &block_grid->interfaces[interfaces_Y_offset + i * (GRID_SIZE_Y-1) + j    ] : NULL;

            struct Block* block = &block_grid->blocks[idx];
            block->state = NewCycle;
            block->cycle = 0;
            block->sweep_num = 0;
            block->must_wait = false;
            block->pos[0] = i;
            block->pos[1] = j;
            block->interfaces[0] = left_int;
            block->interfaces[1] = right_int;
            block->interfaces[2] = bottom_int;
            block->interfaces[3] = top_int;

            // TODO: tmp
//            const int global_pos[] = { block_grid->global_origin[0] + i, block_grid->global_origin[1] + j };
//            for (int s = 0; s < 4; s++) {
//                if (neighbour_ranks[s] == MPI_PROC_NULL) { continue; }
//                enum Side side = s;
//                int axis = s / 2;
//                int axis_offset = (side % 2 == 0) ? -1 : 1;
//                int remote_pos[] = { global_pos[0], global_pos[1] };
//                remote_pos[axis] += axis_offset;
//                if (!(0 <= remote_pos[axis] && remote_pos[axis] < global_grid_size[axis])) {
//                    continue;
//                }
//                init_remote_block(&block_grid->remote_blocks[remote_block_i], comm, side);  // TODO: oops how do we trigger a halo exchange without a pointer to the RemoteBlock?
//                remote_block_i++;
//            }
        }
    }

    int total_assigned_blocks = 0;
    block_grid->num_threads = num_threads;
    block_grid->threads_workload = malloc(sizeof(struct ThreadWorkload) * num_threads);
    for (int tid = 0; tid < num_threads; tid++) {
        init_thread_workload(&block_grid->threads_workload[tid], block_grid->blocks, num_threads, tid);
        total_assigned_blocks += block_grid->threads_workload[tid].num_blocks;
    }

    CIVL_assert(total_assigned_blocks == TOTAL_BLOCKS);

    check_interfaces(block_grid);

    return block_grid;
}

void free_grid(struct BlockGrid* block_grid)
{
    for (int tid = 0; tid < block_grid->num_threads; tid++) {
        free(block_grid->threads_workload[tid].threads_blocks);
    }
    // TODO: tmp
//    for (int i = 0; i < TOTAL_REMOTE_BLOCKS; i++) {
//        free_remote_block(&block_grid->remote_blocks[i]);
//    }
    free(block_grid->threads_workload);
    free(block_grid->blocks);
    free(block_grid->remote_blocks);
    free(block_grid->interfaces);
    free(block_grid);
}

void block_ghost_exchange(struct Block* block)
{
#if DO_HALO_EXCHANGE
    enum Side side;
    enum BlockExchangeState left_exchange_state, right_exchange_state;
    struct BlockInterface *interface, *left_interface, *right_interface;

    side = block->sweep_num % 2 == 0 ? Left : Bottom;
    left_interface = interface = block->interfaces[side];
    if (interface != NULL && !interface->is_done[0]) {
        // block_ghost_exchange between two local blocks
        bool do_xchg = mark_ready_for_exchange(interface, true, &left_exchange_state);
        if (do_xchg) {
            // "do the exchange between the blocks"
            left_exchange_state = exchange_done(interface);
        }
        if (left_exchange_state == XCHG_Done) {
            interface->is_done[0] = true;
        }
    } else {
        left_exchange_state = XCHG_Done;
    }

    side = block->sweep_num % 2 == 0 ? Right : Top;
    right_interface = interface = block->interfaces[side];
    if (interface != NULL && !interface->is_done[1]) {
        // block_ghost_exchange between two local blocks
        bool do_xchg = mark_ready_for_exchange(interface, false, &right_exchange_state);
        if (do_xchg) {
            // "do the exchange between the blocks"
            right_exchange_state = exchange_done(interface);
        }
        if (right_exchange_state == XCHG_Done) {
            interface->is_done[1] = true;
        }
    } else {
        right_exchange_state = XCHG_Done;
    }

    if (left_exchange_state == XCHG_Done && right_exchange_state == XCHG_Done) {
        if (left_interface  != NULL) left_interface ->is_done[0] = false;
        if (right_interface != NULL) right_interface->is_done[1] = false;
        block->must_wait = false;
    } else {
        block->must_wait = true;
    }
#else
    block->must_wait = false;
#endif
}

void block_state_machine(struct Block* block)
{
    bool stop_processing = false;
    for (; !stop_processing;) {
        switch (block->state) {
        case NewCycle: {
            if (block->cycle == global_cycle) {
                block->state = TimeStep;
            } else {
                stop_processing = true;
            }
            break;
        }
        case TimeStep:
        case InitTimeStep: {
            next_time_step(block);
            if (block->must_wait) {
                block->state = InitTimeStep;
                stop_processing = true;
                break;
            } else {
                block->state = NewSweep;
            }
            break;
        }
        case NewSweep: {
            if (block->sweep_num == MAX_SWEEPS-1) {
                block->sweep_num = 0;
                block->state = EndCycle;
            } else {
                block->sweep_num++;
                block->state = EOS;
            }
            break;
        }
        case EOS: { block->state = Exchange; break; }
        case Exchange: {
            block_ghost_exchange(block);
            if (block->must_wait) {
                stop_processing = true;
            } else {
                block->state = Fluxes;
            }
            break;
        }
        case Fluxes:     { block->state = CellUpdate; break; }
        case CellUpdate: { block->state = Remap;      break; }
        case Remap:      { block->state = NewSweep;   break; }
        case EndCycle: {
            block->cycle++;
            block->state = NewCycle;
            stop_processing = true;
            break;
        }
        }
    }
}
