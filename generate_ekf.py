from sympy import *

# states, inputs, and process noise
X=Matrix(symbols('X[0] X[1] X[2] X[3] X[4] X[5] X[6] X[7] X[8] X[9] X[10] X[11] X[12] X[13] X[14] X[15]'))
U=Matrix(symbols('U[0] U[1] U[2] U[3] U[4] U[5]'))
W=Matrix(symbols('W[0] W[1] W[2] W[3] W[4] W[5] W[6] W[7] W[8] W[9] W[10] W[11]'))

# time step
dt=symbols('dt')

# continuous dynamics:
g = 9.81
x,y,z,vx,vy,vz,qw,qx,qy,qz,lx,ly,lz,lp,lq,lr = X
ax,ay,az,p,q,r = U
wx,wy,wz,wp,wq,wr,wbx,wby,wbz,wbp,wbq,wbr = W

quat = Quaternion(qw, qx, qy, qz, norm=1) # does norm 1 automatically renormaize?
a_NED = Quaternion.rotate_point([ax-lx-wx,ay-ly-wy,az-lz-wz], quat)

# https://ahrs.readthedocs.io/en/latest/filters/angular.html#quaternion-derivative
pqr_hat = Matrix([p-lp-wp, q-lq-wq, r-lr-wr])
Omega_pqr = Matrix([
    [0,          -pqr_hat[0], -pqr_hat[1], -pqr_hat[2]],
    [pqr_hat[0],  0,           pqr_hat[2], -pqr_hat[1]],
    [pqr_hat[1], -pqr_hat[2],  0,          +pqr_hat[0]],
    [pqr_hat[2], +pqr_hat[1], -pqr_hat[0], 0],
])
q_dot = 0.5*Omega_pqr * Matrix([quat.a, quat.b, quat.c, quat.d])

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
    #0,0,0,0,0,0
    wbx, wby, wbz, wbp, wbq, wbr
])

# discretized dynamics:
f = X + f_continuous*dt

# output function (measurement model):
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
symmetric_indexing = lambda i,j: int(j*(j+1)/2+i) if j>=i else int(i*(i+1)/2+j)
P = Matrix([[symbols(f'P[{symmetric_indexing(i,j)}]') for j in range(len(X))] for i in range(len(X))])
Q = Matrix([[symbols(f'Q[{i}]') if i==j else '0' for j in range(len(W))] for i in range(len(W))])
R = Matrix([[symbols(f'R[{i}]') if i==j else '0' for j in range(len(h))] for i in range(len(h))])
Z = Matrix([symbols(f'Z[{i}]') for i in range(len(h))])


#%% generate code

# from sympy import symbols, sin
from sympy.codegen.ast import CodeBlock, Assignment

# EKF equations from: https://en.wikipedia.org/wiki/Extended_Kalman_filter: Non-additive noise formulation and equations

# PREDICTION STEP
Xpred = f
Ppred = F*P*F.T + L*Q*L.T
for i in range(6,10):
    Ppred[i, i] += dt*dt*1e-2

# UPDATE STEP
S = H*P*H.T + R
sdim = S.shape[0]
xdim = len(X)

## old way with symbolic inversion
#inv__ = invert_function(sdim)
#iS = inv__(S)
#K = P*H.T*iS

## new way with numerical inversion
HP = H*P
# numerical solution using Cholesky factorization of  S*K^T = HP  takes place here and produces column major K
#                                         columns ,              rows
K = Matrix([[symbols(f'K[{j*xdim + i}]') for j in range(sdim)] for i in range(xdim)])

Xup = X + K*(Z - h)
Pup = (eye(len(X)) - K*H)*P

# assignments
Xpred_assigments = [Assignment(symbols(f'X_new[{i}]'), Xpred[i]) for i in range(len(X))]
Ppred_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Ppred[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated

S_assigments = [Assignment(symbols(f'S[{j*sdim+i}]'), S[i,j]) for i in range(len(Z)) for j in range(len(Z))] # both triangular, full matrix
HP_assigments = [Assignment(symbols(f'HP[{j*sdim+i}]'), HP[i,j]) for j in range(xdim) for i in range(sdim)]

Xup_assigments = [Assignment(symbols(f'X_new[{i}]'), Xup[i]) for i in range(len(X))]
Pup_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Pup[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated


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
s_code = CodeBlock(*S_assigments, *HP_assigments)
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


#%% post-proc code

# replace sin, cos, tan, pow with sinf, cosf, tanf, powf
prediction_code = prediction_code.replace('sin(', 'sinf(')
prediction_code = prediction_code.replace('cos(', 'cosf(')
prediction_code = prediction_code.replace('tan(', 'tanf(')
prediction_code = prediction_code.replace('pow(', 'powf(')
prediction_code = prediction_code.replace('\n', '\n\t')
s_code = s_code.replace('sin(', 'sinf(')
s_code = s_code.replace('cos(', 'cosf(')
s_code = s_code.replace('tan(', 'tanf(')
s_code = s_code.replace('pow(', 'powf(')
s_code = s_code.replace('\n', '\n\t')
update_code = update_code.replace('sin(', 'sinf(')
update_code = update_code.replace('cos(', 'cosf(')
update_code = update_code.replace('tan(', 'tanf(')
update_code = update_code.replace('pow(', 'powf(')
update_code = update_code.replace('\n', '\n\t')


#%% emit code

from os import path, makedirs
from jinja2 import Environment, FileSystemLoader

# set data used in templates
data = {}
data['lenX'] = len(X)
data['lenQ'] = len(W)
data['lenU'] = len(U)
data['lenh'] = len(h)
data['lenTmp'] = max(prediction_code_N_tmps, s_code_N_tmps, update_code_N_tmps)
data['prediction_code'] = prediction_code
data['prepare_gain_code'] = s_code
data['update_code'] = update_code

# make dirs 
home_path = path.dirname(__file__)
code_path = path.join(home_path, 'c_code')
makedirs(code_path, exist_ok=True)

# create files from templates
files = [ "ekf_calc.h", "ekf_calc.c" ]

env = Environment(loader = FileSystemLoader(home_path),
                  trim_blocks=True, lstrip_blocks=True)

for filename in files:
    template = env.get_template(f'{filename}.j2')
    with open(path.join(code_path, filename), 'w') as file:
        file.write(template.render(data))

