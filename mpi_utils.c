#include "mpi.h"

#include "consts.h"

void exchange_halos (int nx_local, double field[nx_local+2][NY][NZ], int left, int right) {
    MPI_Request reqs[4];

    // first, send left slab left
    MPI_Isend(field[1], NY * NZ, MPI_DOUBLE, left, SEND_L, MPI_COMM_WORLD, &reqs[0]);

    // then, send right slab right
    MPI_Isend(field[nx_local], NY * NZ, MPI_DOUBLE, right, SEND_R, MPI_COMM_WORLD, &reqs[1]);

    // then, receive neighboring slabs
    MPI_Irecv(field[nx_local+1], NY * NZ, MPI_DOUBLE, right, SEND_L, MPI_COMM_WORLD, &reqs[2]);
    MPI_Irecv(field[0], NY * NZ, MPI_DOUBLE, left, SEND_R, MPI_COMM_WORLD, &reqs[3]);

    // finally waitall
    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
}
