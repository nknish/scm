#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mpi.h"

#include "consts.h"

void print_slice(int rank, int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    const char* msg);
void exchange_halos(int nx_local, double field[nx_local+2][NY][NZ], int left, int right);
void init_fields(int rank, int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ]);
void compute_advection(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ],
    double dtheta_dt[nx_local+2][NY][NZ]);
void compute_coriolis(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ]);
void apply_tendencies(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ],
    double dtheta_dt[nx_local+2][NY][NZ]);

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
        // halo exchange
        exchange_halos(nx_local, u, left, right);
        exchange_halos(nx_local, v, left, right);
        exchange_halos(nx_local, theta, left, right);

        // advection
        compute_advection(nx_local, u, v, theta, du_dt, dv_dt, dtheta_dt);

        // coriolis
        compute_coriolis(nx_local, u, v, du_dt, dv_dt);

        // numerical time-step
        apply_tendencies(nx_local, u, v, theta, du_dt, dv_dt, dtheta_dt);

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

void init_fields (int rank, int nx_local, double u[nx_local+2][NY][NZ], double v[nx_local+2][NY][NZ], double theta[nx_local+2][NY][NZ]) {
    // initilize u with a special pattern
    for (int i = 0; i < nx_local+2; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                // Interior: use rank*3 + local i for visibility
                if (i == 0 || i == nx_local+1) {
                    u[i][j][k] = -1.0; // initial halo marker
                } else {
                    u[i][j][k] = rank*3 + i;
                }
            }
        }
    }

    // initialize v and theta uniformly
    for (int i = 0; i < nx_local+2; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                v[i][j][k] = (double) rank;
                theta[i][j][k] = (double) rank;
            }
        }
    }
}

void compute_advection(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ],
    double dtheta_dt[nx_local+2][NY][NZ]) {

    // advection: du/dt + u*(du/dx) + v(dv/dy) + w(dw/dz) = 0, ignoring w (vertical)
    for (int i = 1; i < nx_local+1; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                // periodic boundaries at poles
                // consider instead averaging, e.g. the north pole
                // (above each max-lat grid box) is the average of all
                // max-lat grid boxes
                int jp = (j + 1) % NY;
                int jm = (j - 1 + NY) % NY;

                // calculate du/dt via advection pde
                double du_dx = (u[i+1][j][k] - u[i-1][j][k]) / (2 * DX);
                double du_dy = (u[i][jp][k] - u[i][jm][k]) / (2 * DY);
                du_dt[i][j][k] = - (u[i][j][k]*du_dx + v[i][j][k]*du_dy); 

                // same for dv/dt
                double dv_dx = (v[i+1][j][k] - v[i-1][j][k]) / (2 * DX);
                double dv_dy = (v[i][jp][k] - v[i][jm][k]) / (2 * DY);
                dv_dt[i][j][k] = - (u[i][j][k]*dv_dx + v[i][j][k]*dv_dy); 

                // and for dtheta/dt
                double dtheta_dx = (theta[i+1][j][k] - theta[i-1][j][k]) / (2 * DX);
                double dtheta_dy = (theta[i][jp][k] - theta[i][jm][k]) / (2 * DY);
                dtheta_dt[i][j][k] = - (u[i][j][k]*dtheta_dx + v[i][j][k]*dtheta_dy); 
            }
        }
    }
}

void compute_coriolis(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ]) {

    // coriolis: subtract f*v from du/dt, add f*u to dv/dt
    double f = 1e-4;
    for (int i = 1; i < nx_local+1; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                du_dt[i][j][k] -= f * v[i][j][k];
                dv_dt[i][j][k] += f * u[i][j][k];
            }
        }
    }
}

void apply_tendencies(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ],
    double dtheta_dt[nx_local+2][NY][NZ]) {

    // euler's method (for now)
    for (int i = 1; i < nx_local+1; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                u[i][j][k] += DT * du_dt[i][j][k];
                v[i][j][k] += DT * dv_dt[i][j][k];
                theta[i][j][k] += DT * dtheta_dt[i][j][k];
            }
        }
    }
}
