#ifndef EULER_H
#define EULER_H

void rk2_step(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ],
    double k1_u[nx_local+2][NY][NZ],
    double k1_v[nx_local+2][NY][NZ],
    double k1_theta[nx_local+2][NY][NZ],
    double k2_u[nx_local+2][NY][NZ],
    double k2_v[nx_local+2][NY][NZ],
    double k2_theta[nx_local+2][NY][NZ],
    double u_tmp[nx_local+2][NY][NZ],
    double v_tmp[nx_local+2][NY][NZ],
    double theta_tmp[nx_local+2][NY][NZ],
    int left, int right);

#endif
