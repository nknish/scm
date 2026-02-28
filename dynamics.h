#ifndef DYNAMICS_H
#define DYNAMICS_H

void init_fields (int rank, int nx_local,
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

#endif
