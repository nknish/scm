#include "consts.h"
#include "dynamics.h"

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
    int left, int right) {

    compute_rhs(nx_local, u, v, theta, k1_u, k1_v, k1_theta, left, right);
    for (int i = 1; i < nx_local+1; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                u_tmp[i][j][k] = u[i][j][k] + DT * k1_u[i][j][k];
                v_tmp[i][j][k] = v[i][j][k] + DT * k1_v[i][j][k];
                theta_tmp[i][j][k] = theta[i][j][k] + DT * k1_theta[i][j][k];
            }
        }
    }
    compute_rhs(nx_local, u_tmp, v_tmp, theta_tmp, k2_u, k2_v, k2_theta, left, right);
    for (int i = 1; i < nx_local+1; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                u[i][j][k] += DT/2 * (k1_u[i][j][k] + k2_u[i][j][k]);
                v[i][j][k] += DT/2 * (k1_v[i][j][k] + k2_v[i][j][k]);
                theta[i][j][k] += DT/2 * (k1_theta[i][j][k] + k2_theta[i][j][k]);
            }
        }
    }
}


