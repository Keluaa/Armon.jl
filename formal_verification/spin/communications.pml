
#ifndef USE_MPI
#define USE_MPI           0
#endif

#ifndef MPI_CHANNEL_SIZE
#define MPI_CHANNEL_SIZE  10
#endif

// MPI neighbours
chan subdomain_neighbours_dt_reduction     = [MPI_CHANNEL_SIZE] of { byte };
chan subdomain_neighbours_halo_exchange[4] = [MPI_CHANNEL_SIZE] of { byte };

// TODO
