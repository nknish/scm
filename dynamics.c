#include <math.h>

#include "consts.h"
#include "mpi_utils.h"

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

void zero_tendencies(int nx_local,
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ],
    double dtheta_dt[nx_local+2][NY][NZ]) {

    for (int i = 0; i < nx_local + 2; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                du_dt[i][j][k] = 0;
                dv_dt[i][j][k] = 0;
                dtheta_dt[i][j][k] = 0;
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
    for (int i = 1; i < nx_local+1; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                du_dt[i][j][k] -= F_CORIOLIS * v[i][j][k];
                dv_dt[i][j][k] += F_CORIOLIS * u[i][j][k];
            }
        }
    }
}

void compute_rhs(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ],
    double dtheta_dt[nx_local+2][NY][NZ],
    int left,
    int right) {

    zero_tendencies(nx_local, du_dt, dv_dt, dtheta_dt);

    // comunication
    exchange_halos(nx_local, u, left, right);
    exchange_halos(nx_local, v, left, right);
    exchange_halos(nx_local, theta, left, right);

    // computation
    compute_advection(nx_local, u, v, theta, du_dt, dv_dt, dtheta_dt);
    compute_coriolis(nx_local, u, v, du_dt, dv_dt);
}
