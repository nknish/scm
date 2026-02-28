#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"

void print_slice (int rank, int nx_local, double u[nx_local+2][NY][NZ], double v[nx_local+2][NY][NZ], const char* msg) {
    //arbitrary point on the grid
    int x = 12;
    int y = 15;
    int z = 2;
    int real_x = x + nx_local * rank;
    double mag = sqrt(pow(u[x][y][z], 2) + pow(v[x][y][z], 2));
    printf("Rank %d at(%d,%d,%d): %s\t", rank, real_x, y, z, msg);
    printf("u: %6.1f \tv: %6.1f \t|u|: %6.1f", u[x][y][z], v[x][y][z], mag);
    printf("\n");
    fflush(stdout);
}

