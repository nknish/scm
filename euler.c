#include "consts.h"

void time_integrate(int nx_local,
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
