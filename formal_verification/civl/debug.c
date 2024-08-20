
#include <stdio.h>

#include "debug.h"


const char* side_to_str(enum Side side)
{
    switch (side) {
    case Left:   return "Left";
    case Right:  return "Right";
    case Bottom: return "Bottom";
    case Top:    return "Top";
    }
}

const char* block_state_to_str(enum BlockState state)
{
    switch (state) {
    case NewCycle:     return "NewCycle";
    case TimeStep:     return "TimeStep";
    case InitTimeStep: return "InitTimeStep";
    case NewSweep:     return "NewSweep";
    case EOS:          return "EOS";
    case Exchange:     return "Exchange";
    case Fluxes:       return "Fluxes";
    case CellUpdate:   return "CellUpdate";
    case Remap:        return "Remap";
    case EndCycle:     return "EndCycle";
    }
}

const char* time_step_state_to_str(enum TimeStepState time_step_state)
{
    switch (time_step_state) {
    case DT_Ready:          return "Ready";
    case DT_AllContributed: return "AllContributed";
    case DT_DoingMPI:       return "DoingMPI";
    case DT_WaitingForMPI:  return "WaitingForMPI";
    case DT_Done:           return "Done";
    }
}

const char* block_xchg_state_to_str(enum BlockExchangeState block_exchange_state)
{
    switch (block_exchange_state) {
    case XCHG_NotReady:   return "NotReady";
    case XCHG_InProgress: return "InProgress";
    case XCHG_Done:       return "Done";
    }
}

void print_interface(struct BlockInterface* interface)
{
    if (interface == NULL) {
        printf("none");
    } else {
        enum BlockExchangeState bint_state;
        byte flags;
        block_interface_state(interface, &bint_state, &flags);
        printf("addr=%p, left=%d, right=%d, state=%s, flags=%d",
               interface, interface->is_done[0], interface->is_done[1], block_xchg_state_to_str(bint_state), flags);
    }
}

void print_block(struct Block* block)
{
    printf("block at (%d,%d), state=%s, cycle=%d, sweep=%d, int={",
           block->pos[0], block->pos[1], block_state_to_str(block->state), block->cycle, block->sweep_num);
    for (int side = 0; side < 4; side++) {
        printf("\n    - [%s]={", side_to_str(side));
        print_interface(block->interfaces[side]);
        printf("}");
    }
    printf("\n   }");
}

void print_block_grid(struct BlockGrid* block_grid)
{
    printf("Block grid, (%d x %d), %d threads:\n", GRID_SIZE_X, GRID_SIZE_Y, block_grid->num_threads);
    for (int block_i = 0; block_i < TOTAL_BLOCKS; block_i++) {
        printf("  ");
        print_block(block_grid->blocks + block_i);
        printf("\n");
    }
}

void print_global_state()
{
    printf("Global state:\n");
    printf(" - global_dt_state  = %s\n", time_step_state_to_str(global_dt_state));
    printf(" - dt_contributions = %d\n", dt_contributions);
    printf(" - global_cycle     = %d\n", global_cycle);
}

void pretty_print_solver_state(struct BlockGrid* block_grid)
{
    printf(" === SOLVER STATE ===\n");
    print_block_grid(block_grid);
    printf("\n");
    print_global_state();
}

#ifdef _CIVL

bool check_interfaces(struct BlockGrid* block_grid) { return true; }

#else

static int find_interface(struct BlockGrid* block_grid, struct BlockInterface* interface)
{
    if (interface == NULL) return -1;
    for (int i = 0; i < TOTAL_INTERFACES; i++) {
        if (&block_grid->interfaces[i] == interface) return i;
    }
    CIVL_assert(false);
    return -2;
}

bool check_interfaces(struct BlockGrid* block_grid)
{
    struct Block* connected_blocks[TOTAL_INTERFACES][2];
    for (int i = 0; i < TOTAL_INTERFACES; i++) {
        connected_blocks[i][0] = connected_blocks[i][1] = NULL;
    }

    bool ok = true;
    for (int tid = 0; tid < block_grid->num_threads; tid++) {
        struct ThreadWorkload* thread_workload = &block_grid->threads_workload[tid];
        for (int block_i = 0; block_i < thread_workload->num_blocks; block_i++) {
            struct Block* block = thread_workload->threads_blocks[block_i];
            for (int side = 0; side < 4; side++) {
                int int_idx = find_interface(block_grid, block->interfaces[side]);
                if (int_idx < 0) continue;
                if (connected_blocks[int_idx][side % 2] == NULL) {
                    connected_blocks[int_idx][side % 2] = block;
                } else {
                    ok = false;
                    printf("Interface %d already connected to (%d,%d) along side %d\n",
                           int_idx, connected_blocks[int_idx][side % 2]->pos[0], connected_blocks[int_idx][side % 2]->pos[1], side);
                }
            }
        }
    }

    for (int i = 0; i < TOTAL_INTERFACES; i++) {
        bool int_ok = true;
        for (int s = 0; s < 2; s++) {
            if (connected_blocks[i][s] == NULL) {
                printf("Interface %d not connected along side %d\n", i, s);
                ok = false;
                int_ok = false;
                continue;
            }
        }
        if (!int_ok) continue;

        int left_pos[2]  = { connected_blocks[i][0]->pos[0], connected_blocks[i][0]->pos[1] };
        int right_pos[2] = { connected_blocks[i][1]->pos[0], connected_blocks[i][1]->pos[1] };
        printf("Interface %d connected to (%d,%d) and (%d,%d)\n", i, left_pos[0], left_pos[1], right_pos[0], right_pos[1]);
    }

    return ok;
}

#endif //_CIVL
