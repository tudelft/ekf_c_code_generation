from sympy import *

# states, inputs, and process noise
X=Matrix(symbols('X[0] X[1] X[2] X[3] X[4] X[5] X[6] X[7] X[8] X[9] X[10] X[11] X[12] X[13] X[14] X[15]'))
U=Matrix(symbols('U[0] U[1] U[2] U[3] U[4] U[5]'))
W=Matrix(symbols('W[0] W[1] W[2] W[3] W[4] W[5]'))

# time step
dt=symbols('dt')

# continuous dynamics:
g = 9.81
x,y,z,vx,vy,vz,qw,qx,qy,qz,lx,ly,lz,lp,lq,lr = X
ax,ay,az,p,q,r = U
wx,wy,wz,wp,wq,wr = W

q = Quaternion(qw, qx, qy, qz, norm=1) # does norm 1 automatically renormaize?
a_NED = Quaternion.rotate_point([ax-lx-wx,ay-ly-wy,az-lz-wz], q)

# https://ahrs.readthedocs.io/en/latest/filters/angular.html#quaternion-derivative
pqr_hat = Matrix([p-lp-wp, q-lq-wq, q-lq-wq])
Omega_pqr = Matrix([
    [0,          -pqr_hat[0], -pqr_hat[1], -pqr_hat[2]],
    [pqr_hat[0],  0,           pqr_hat[2], -pqr_hat[1]],
    [pqr_hat[1], -pqr_hat[2],  0,          +pqr_hat[0]],
    [pqr_hat[2], +pqr_hat[1], -pqr_hat[0], 0],
])
q_dot = 0.5*Omega_pqr * Matrix([q.a, q.b, q.c, q.d])

f_continuous = Matrix([
    vx,
    vy,
    vz,
    a_NED[0],
    a_NED[1],
    a_NED[2] + g,
    q_dot[0],
    q_dot[1],
    q_dot[2],
    q_dot[3],
    0,0,0,0,0,0
])


# discretized dynamics:
f = X + f_continuous*dt

# output function (measurement model):
#use_phi, use_theta, use_psi = symbols('ekf_use_phi ekf_use_theta ekf_use_psi') # 0 or 1
h = Matrix([x,y,z,qw,qx,qy,qz])

# matrices:
F = f.jacobian(X)
L = f.jacobian(W)
H = h.jacobian(X)

# substitute W with 0
f = f.subs([(w,0) for w in W])
F = F.subs([(w,0) for w in W])
L = L.subs([(w,0) for w in W])
H = H.subs([(w,0) for w in W])

# extra matrices
symmetric_indexing = lambda i,j: int(i*(i+1)/2+j) if i>=j else int(j*(j+1)/2+i)
P = Matrix([[symbols(f'P[{symmetric_indexing(i,j)}]') for j in range(len(X))] for i in range(len(X))])
Q = Matrix([[symbols(f'Q[{i}]') if i==j else '0' for j in range(len(W))] for i in range(len(W))])
R = Matrix([[symbols(f'R[{i}]') if i==j else '0' for j in range(len(h))] for i in range(len(h))])
Z = Matrix([symbols(f'Z[{i}]') for i in range(len(h))])


#%%
# function for inverting matrix with symbolic elements that is quicker
def invert_function(num):
    M = Matrix([[symbols(f'M[{i}{j}]', positive=(i==j)) for j in range(num)] for i in range(num)])
    M_inv = M.inverse_CH()
    return lambda M_: M_inv.subs({M[i,j]:M_[i,j] for i in range(num) for j in range(num)})



#%%

# from sympy import symbols, sin
from sympy.codegen.ast import CodeBlock, Assignment

# EKF equations from: https://en.wikipedia.org/wiki/Extended_Kalman_filter: Non-additive noise formulation and equations

# PREDICTION STEP
Xpred = f
Ppred = F*P*F.T + L*Q*L.T

# UPDATE STEP
S = H*P*H.T + R
sdim = S.shape[0]
xdim = len(X)

## old way with symbolic inversion
#inv__ = invert_function(sdim)
#iS = inv__(S)
#K = P*H.T*iS

## new way with numerical inversion
PHT = P*H.T
# numerical solution using Cholesky factorization of  K*S = PHT  takes place here and produces K
K = Matrix([[symbols(f'K[{i}{j}]') for j in range(sdim)] for i in range(xdim)])

Xup = X + K*(Z - h)
Pup = (eye(len(X)) - K*H)*P

# assignments
Xpred_assigments = [Assignment(symbols(f'X_new[{i}]'), Xpred[i]) for i in range(len(X))]
Ppred_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Ppred[i,j]) for i in range(len(X)) for j in range(len(X)) if i >= j] # only the lower diagonal will be calculated

S_assigments = [Assignment(symbols(f'S[{symmetric_indexing(i,j)}]'), S[i,j]) for i in range(len(Z)) for j in range(len(Z)) if i >= j] # only the lower diagonal will be calculated
PHT_assigments = [Assignment(symbols(f'PHT[{j*xdim+i}]'), PHT[i,j]) for i in range(xdim) for j in range(sdim)]

Xup_assigments = [Assignment(symbols(f'X_new[{i}]'), Xup[i]) for i in range(len(X))]
Pup_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Pup[i,j]) for i in range(len(X)) for j in range(len(X)) if i >= j] # only the lower diagonal will be calculated

# PREDICTION STEP
print('PREDICTION:')
# code generation
prediction_code = CodeBlock(*Xpred_assigments, *Ppred_assigments)
# common subexpression elimination with tmp variables
print('CSE')
prediction_code = prediction_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
# simplify
print('SIMPLIFY')
prediction_code = prediction_code.simplify()
# count number of tmp variables
prediction_code_N_tmps = len([s for s in prediction_code.left_hand_sides if s.name.startswith('tmp')])
# generate C code
print('CCODE')
prediction_code = ccode(prediction_code)


# Prepare numerical solution for K
print('S and PHT MATRICES FOR SOLVING K:')
# code generation
s_code = CodeBlock(*S_assigments, *PHT_assigments)
# common subexpression elimination with tmp variables
print('CSE')
s_code = s_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
# simplify
print('SIMPLIFY')
s_code = s_code.simplify()
# count number of tmp variables
s_code_N_tmps = len([s for s in s_code.left_hand_sides if s.name.startswith('tmp')])
# generate C code
print('CCODE')
s_code = ccode(s_code)


# UPDATE STEP
print('UPDATE:')
# code generation
update_code = CodeBlock(*Xup_assigments, *Pup_assigments)
# common subexpression elimination with tmp variables
print('CSE')
update_code = update_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
# simplify
print('SIMPLIFY')
update_code = update_code.simplify()
# count number of tmp variables
update_code_N_tmps = len([s for s in update_code.left_hand_sides if s.name.startswith('tmp')])
# generate C code
print('CCODE')
update_code = ccode(update_code)



#%% 

# replace sin, cos, tan, pow with sinf, cosf, tanf, powf
prediction_code = prediction_code.replace('sin(', 'sinf(')
prediction_code = prediction_code.replace('cos(', 'cosf(')
prediction_code = prediction_code.replace('tan(', 'tanf(')
prediction_code = prediction_code.replace('pow(', 'powf(')
update_code = update_code.replace('sin(', 'sinf(')
update_code = update_code.replace('cos(', 'cosf(')
update_code = update_code.replace('tan(', 'tanf(')
update_code = update_code.replace('pow(', 'powf(')

print('Prediction code:')
print(prediction_code)
print('Update code:')
print(update_code)


#%%

# create c_code folder (if it doesn't exist)
import os
if not os.path.exists('c_code'):
    os.makedirs('c_code')
    
header = f'''
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This file is automatically generated in "Generate Kalman Filter C Code.ipynb" from https://github.com/tudelft/kalman_filter   //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EKF_CALC_H
#define EKF_CALC_H

#include <stdint.h>

#define N_STATES {len(X)}
#define N_INPUTS {len(U)}
#define N_MEASUREMENTS {len(h)}

// set to 1 to use phi, theta, psi measurements
float ekf_use_phi;
float ekf_use_theta;
float ekf_use_psi;

// getters
float* ekf_get_X();     // get state vector
float* ekf_get_P();     // get covariance matrix (lower diagonal)

#define ekf_P_index(i,j) ((i>=j) ? ekf_get_P()[i*(i+1)/2+j] : ekf_get_P()[j*(j+1)/2+i])
#define ekf_X_index(i) ekf_get_X()[i]

// setters
void ekf_set_Q(float Q[N_INPUTS]);                    // set process noise covariance matrix diagonal
void ekf_set_R(float R[N_MEASUREMENTS]);              // set measurement noise covariance matrix diagonal
void ekf_set_X(float X0[N_STATES]);                   // set state vector
void ekf_set_P_diag(float P_diag[N_STATES]);          // set covariance matrix diagonal
void ekf_set_P(float P0[N_STATES*(N_STATES+1)/2]);    // set covariance matrix (lower diagonal)

// prediction and update functions
void ekf_predict(float U[N_INPUTS], float dt);
void ekf_update(float Z[N_MEASUREMENTS]);


#endif // EKF_CALC_H
'''

prediction_code_ = prediction_code.replace('\n','\n\t')
update_code_ = update_code.replace('\n','\n\t')

source = f'''
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This file is automatically generated in "Generate Kalman Filter C Code.ipynb" from https://github.com/tudelft/kalman_filter   //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "ekf_calc.h"
#include <math.h>

void pprz_cholesky_float(num_t **out, num_t **in, num_t* inv_diag, int n)
{
  int i, j, k;

  for (i = 0; i < n; i++) {
    for (j = 0; j < (i + 1); j++) {
      num_t s = 0;
      for (k = 0; k < j; k++) {
        s += out[i][k] * out[j][k];
      }
      if (i == j) {
        out[i][j] = sqrtf(in[i][i] - s);
        inv_diag[j] = (out[i][j] != 0) ? (1.0 / out[i][j]) : 0.0;
      } else {
        out[i][j] = inv_diag[j] * (in[i][j] - s);
      }
    }
  }
}

void cholesky_solve(num_t **L, num_t* inv_diag, int n, num_t *b, num_t *x) {
    // Antoine Drouin, 2007, modified
	int j,k;
	num_t t;

    for(j = 0 ; j < n ; j++) { // solve Ly=b
        t = b[j];
        for(k = j - 1 ; k >= 0 ; k--)
            t -= L[j][k] * x[k];
        x[j] = t*inv_diag[j];
    }
    for(j = n - 1 ; j >= 0 ; j--) { // solve Ltx=y
        t = x[j];
        for(k = j + 1 ; k < n ; k++)
            t -= L[k][j] * x[k];
        x[j] = t*inv_diag[j];
    }
}

float ekf_use_phi = 1;
float ekf_use_theta = 1;
float ekf_use_psi = 1;

float ekf_Q[N_INPUTS];         // Kalman filter process noise covariance matrix (diagonal)
float ekf_R[N_MEASUREMENTS];   // Kalman filter measurement noise covariance matrix (diagonal)

// state
float ekf_X[N_STATES];
float ekf_X_new[N_STATES];

// covariance matrix (lower diagonal) P[i,j] = P_lower_diagonal[i*(i+1)/2+j] (if i>=j)
float ekf_P_lower_diagonal[N_STATES*(N_STATES+1)/2];
float ekf_P_lower_diagonal_new[N_STATES*(N_STATES+1)/2];

// temporary variables
float tmp[{max(prediction_code_N_tmps,update_code_N_tmps)}];

// pointers
float *X = ekf_X;
float *X_new = ekf_X_new;

float *P = ekf_P_lower_diagonal;
float *P_new = ekf_P_lower_diagonal_new;

// renaming
float *Q = ekf_Q;
float *R = ekf_R;

// pointer for swapping
float *swap_ptr;

float* ekf_get_X() {{
    return X;
}}

float* ekf_get_P() {{
    return P;
}}

void ekf_set_Q(float Q[N_INPUTS]) {{
    for (int i=0; i<N_INPUTS; i++) {{
        ekf_Q[i] = Q[i];
    }}
}}

void ekf_set_R(float R[N_MEASUREMENTS]) {{
    for (int i=0; i<N_MEASUREMENTS; i++) {{
        ekf_R[i] = R[i];
    }}
}}

void ekf_set_X(float X0[N_STATES]) {{
    for (int i=0; i<N_STATES; i++) {{
        X[i] = X0[i];
    }}
}}

void ekf_set_P_diag(float P_diag[N_STATES]) {{
    // set P to zeros
    for (int i=0; i<N_STATES*(N_STATES+1)/2; i++) {{
        P[i] = 0.;
    }}
    // set diagonal
    for (int i=0; i<N_STATES; i++) {{
        P[i*(i+1)/2+i] = P_diag[i];
    }}
}}

void ekf_set_P(float P0[N_STATES*(N_STATES+1)/2]) {{
    for (int i=0; i<N_STATES*(N_STATES+1)/2; i++) {{
        P[i] = P0[i];
    }}
}}

void ekf_predict(float U[N_INPUTS], float dt) {{
    // PREDICTION STEP X_new, P_new = ...
    {prediction_code_}

    // swap X, X_new and P, P_new pointers
    swap_ptr = X;
    X = X_new;
    X_new = swap_ptr;

    swap_ptr = P;
    P = P_new;
    P_new = swap_ptr;
}}

void ekf_update(float Z[N_MEASUREMENTS]) {{
    // UPDATE STEP X_new, P_new = ...
    {update_code_}

    // swap X, X_new and P, P_new pointers
    swap_ptr = X;
    X = X_new;
    X_new = swap_ptr;

    swap_ptr = P;
    P = P_new;
    P_new = swap_ptr;
}}
'''

with open('c_code/ekf_calc.h','w') as file:
    file.write(header)
with open('c_code/ekf_calc.c','w') as file:
    file.write(source)

