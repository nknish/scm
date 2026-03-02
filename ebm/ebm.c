#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define S 1361.0
#define ALPHA 0.4
#define SIGMA 5.670374419e-8
#define C 1.0e7

double calc_dT_dt(double T);

int main (int argc, char* argv[]) {
    double T = 0.0;
    double dt = 60.0;  // 1-minute timestep
    double dT = calc_dT_dt(T);
    for (int t = 0; fabs(dT*dt) > 0.000001; t++) {
        if (t % 3600 == 0) {
            printf("after %d days:   ", t / 3600);
            printf("temp: %f,   ", T);
            printf("dT: %f\n", dT);
        }
        T += dT * dt;
        dT = calc_dT_dt(T);
    }

    printf("final temp: %f\n", T);
    return 0;
}

double calc_dT_dt(double T) {
    double absorbed = S * (1 - ALPHA) / 4;
    double emitted = SIGMA * pow(T, 4);
    return 1 / C * (absorbed - emitted);
}
