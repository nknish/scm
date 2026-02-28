#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

#include "consts.h"
#include "dynamics.h"
#include "io.h"

int main (int argc, char* argv[]) {

    int rank;
    int num_ranks;
    int nx_local;
    double t_start, t_end, t_elapsed;

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    if (NX % num_ranks != 0) {
        if (rank == 0) {
            fprintf(stderr, "aborting. num_ranks must be a factor of %d.\n", NX);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // handle dimensions
    nx_local = NX / num_ranks;
    // int local_dim = nx_local * NY * NZ;

    // find neighbors (only in x direction, periodic)
    int left = ((rank - 1) + num_ranks) % num_ranks;
    int right = (rank + 1) % num_ranks;

    // allocate data arrays (+2 in x dimension for halo zones) and tendency (d_field/dt) array
    double (*u)[NY][NZ] = malloc(sizeof(double[nx_local+2][NY][NZ]));
    double (*v)[NY][NZ] = malloc(sizeof(double[nx_local+2][NY][NZ]));
    double (*theta)[NY][NZ] = malloc(sizeof(double[nx_local+2][NY][NZ]));
    double (*du_dt)[NY][NZ] = malloc(sizeof(double[nx_local+2][NY][NZ]));
    double (*dv_dt)[NY][NZ] = malloc(sizeof(double[nx_local+2][NY][NZ]));
    double (*dtheta_dt)[NY][NZ] = malloc(sizeof(double[nx_local+2][NY][NZ]));
    if (!u || !v || !theta || !du_dt || !dv_dt || !dtheta_dt) {
        fprintf(stderr, "memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // initialize fields
    init_fields(rank, nx_local, u, v, theta);
    print_slice(rank, nx_local, u, v, "initial conditions, before first exchange");

    // sync up and start timer
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime();

    for (int t = 0; t < NUM_T_STEPS; t++) {
        // new version
        compute_rhs(nx_local, u, v, theta, du_dt, dv_dt, dtheta_dt, left, right);
        time_integrate(nx_local, u, v, theta, du_dt, dv_dt, dtheta_dt);

        // sanity check
        if (t % 200 == 0 && rank == 0)
            print_slice(rank, nx_local, u, v, "in loop");
    }

    // wrap up timer and output
    MPI_Barrier(MPI_COMM_WORLD);
    t_end = MPI_Wtime();
    t_elapsed = t_end - t_start;

    print_slice(rank, nx_local, u, v, "final state");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
        printf("total time for %d timesteps: %f\n", NUM_T_STEPS, t_elapsed);

    // free memory and finalize
    free(u);
    free(v);
    free(theta);
    free(du_dt);
    free(dv_dt);
    free(dtheta_dt);

    MPI_Finalize();
    return 0;
}
