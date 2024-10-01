#ifndef ARMON_CIVL_MPI_COMMS_H
#define ARMON_CIVL_MPI_COMMS_H

#include <mpi.h>

#include "blocks.h"

struct CommunicationData;

void init_remote_block(struct RemoteBlock* remote_block, MPI_Comm comm, enum Side side);
void free_remote_block(struct RemoteBlock* remote_block);
void init_exchange(struct CommunicationData* comm_data, MPI_Comm comm, int rank, int tag);

void mpi_check(int error_code, const char* file, int line, const char* call);
#define MPI_CHECK(expr) mpi_check((expr), __FILE_NAME__, __LINE__, #expr)

#endif //ARMON_CIVL_MPI_COMMS_H
