
#include <stdlib.h>

#include "mpi_comms.h"

void mpi_check(int error_code, const char* file, int line, const char* call)
{
#ifndef _CIVL
    if (error_code == MPI_SUCCESS) { return; }
    fprintf(stderr, "MPI call at %s line %d returned %d: %s\n", file, line, error_code, call);
    MPI_Abort(MPI_COMM_WORLD, error_code);
#else
    CIVL_assert(error_code == MPI_SUCCESS);
#endif
}


struct CommunicationData {
    int buffer_size;
    double* send_buffer;
    double* recv_buffer;
    MPI_Request send_request;
    MPI_Request recv_request;
};


void init_remote_block(struct RemoteBlock* remote_block, MPI_Comm comm, enum Side side)
{
    int axis = ((int) side) / 2;  // 0: X, 1: Y
    int next_axis = (axis + 1) % 2;
    int tag = remote_block->pos[next_axis];  // unique along blocks of 'side', and common with the remote process

    int buffer_size = 50000;
    remote_block->comm_data->send_buffer = malloc(sizeof(double) * buffer_size);
    remote_block->comm_data->recv_buffer = malloc(sizeof(double) * buffer_size);

    init_exchange(remote_block->comm_data, comm, remote_block->rank, tag);
}

void free_remote_block(struct RemoteBlock* remote_block)
{
    free(remote_block->comm_data->send_buffer);
    free(remote_block->comm_data->recv_buffer);
}

void init_exchange(struct CommunicationData* comm_data, MPI_Comm comm, int rank, int tag)
{
    MPI_Send_init(
        comm_data->send_buffer, comm_data->buffer_size,
        MPI_DOUBLE, rank, tag, comm,
        &comm_data->send_request
    );

    MPI_Recv_init(
        comm_data->recv_buffer, comm_data->buffer_size,
        MPI_DOUBLE, rank, tag, comm,
        &comm_data->recv_request
    );
}
