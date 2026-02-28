#ifndef EULER_H
#define EULER_H

void apply_tendencies(int nx_local,
    double u[nx_local+2][NY][NZ],
    double v[nx_local+2][NY][NZ],
    double theta[nx_local+2][NY][NZ],
    double du_dt[nx_local+2][NY][NZ],
    double dv_dt[nx_local+2][NY][NZ],
    double dtheta_dt[nx_local+2][NY][NZ]);

#endif
