#ifndef ARMON_CIVL_BLOCKS_H
#define ARMON_CIVL_BLOCKS_H

#include <mpi.h>
#include "vars.h"
#include "utils.h"
#include "block_interface.h"


enum Side {
    Left   = 0,
    Right  = 1,
    Bottom = 2,
    Top    = 3
};


enum BlockState {
    NewCycle,
    TimeStep,
    InitTimeStep,
    NewSweep,
    EOS,
    Exchange,
    Fluxes,
    CellUpdate,
    Remap,
    EndCycle
};


struct Block {
    enum BlockState state;
    byte cycle;  // local cycle of the block
    byte sweep_num;
    bool must_wait;  // if the block is waiting for other blocks
    byte pos[2];
    struct BlockInterface* interfaces[4];
};


struct CommunicationData;

struct RemoteBlock {
    int rank;
    struct CommunicationData* comm_data;
    byte pos[2];
    struct BlockInterface* interface;
};


struct BlockGrid {
    byte global_origin[2];
    struct Block* blocks;
    struct RemoteBlock* remote_blocks;
    struct BlockInterface* interfaces;
    int num_threads;
    struct ThreadWorkload {
        int num_blocks;
        struct Block** threads_blocks;
    } *threads_workload;
};


enum TimeStepState {
    DT_Ready,
    DT_AllContributed,
    DT_DoingMPI,
    DT_WaitingForMPI,
    DT_Done
};

// Global state
extern AtomicVar enum TimeStepState global_dt_state;
extern AtomicVar byte dt_contributions;  // number of blocks which contributed to the time step calculation for this cycle
extern byte global_cycle;  // the current cycle of the whole solver


struct BlockGrid* init_grid(MPI_Comm comm, const int* rank_pos, const int* neighbour_ranks, int num_threads);
void free_grid(struct BlockGrid* block_grid);
void block_ghost_exchange(struct Block* block);
void block_state_machine(struct Block* block);

#endif //ARMON_CIVL_BLOCKS_H
