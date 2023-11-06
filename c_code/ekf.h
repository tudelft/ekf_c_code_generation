
#ifndef EKF_H
#define EKF_H

#include <stdint.h>

#define N_STATES 15
#define N_INPUTS 6
#define N_MEASUREMENTS 4

void ekf_init(float Q_diag[N_INPUTS], float R_diag[N_MEASUREMENTS], float X0[N_STATES], float P_diag0[N_STATES]);
float* ekf_get_state();
void ekf_predict(float U[N_INPUTS], float dt);
void ekf_update(float Z[N_MEASUREMENTS]);

#endif
