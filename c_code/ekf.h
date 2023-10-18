
#ifndef EKF_H
#define EKF_H

#include <stdint.h>

#define N_STATES 15
#define N_INPUTS 6
#define N_MEASUREMENTS 6

extern float X[N_STATES];
extern float P[N_STATES][N_STATES];
extern float Q_diag[N_INPUTS];
extern float R_diag[N_MEASUREMENTS];

void ekf_predict(float U[N_INPUTS], float dt);
void ekf_update(float Z[N_MEASUREMENTS]);

#endif
