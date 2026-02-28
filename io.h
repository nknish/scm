#ifndef IO_H
#define IO_H

void print_slice (int rank, int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    const char* msg);

#endif
