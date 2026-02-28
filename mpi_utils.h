#ifndef MPI_UTILS_H
#define MPI_UTILS_H

void exchange_halos (int nx_local, double field[nx_local+2][NY][NZ], int left, int right);

#endif
