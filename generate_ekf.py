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

quat = Quaternion(qw, qx, qy, qz) # does norm 1 automatically renormaize?
a_NED = Quaternion.rotate_point([ax-lx-wx,ay-ly-wy,az-lz-wz], quat)

# https://ahrs.readthedocs.io/en/latest/filters/angular.html#quaternion-derivative
#pqr_hat = Matrix([p-lp-wp, q-lq-wq, r-lr-wr])
#Omega_pqr = Matrix([
#    [0,          -pqr_hat[0], -pqr_hat[1], -pqr_hat[2]],
#    [pqr_hat[0],  0,           pqr_hat[2], -pqr_hat[1]],
#    [pqr_hat[1], -pqr_hat[2],  0,          +pqr_hat[0]],
#    [pqr_hat[2], +pqr_hat[1], -pqr_hat[0], 0],
#])
#q_dot = 0.5*Omega_pqr * Matrix([quat.a, quat.b, quat.c, quat.d])

pqr_hat = Quaternion(0, p-lp-wp, q-lq-wq, r-lr-wr)
q_dot = 0.5 * quat * pqr_hat

f_continuous = Matrix([
    vx,
    vy,
    vz,
    a_NED[0],
    a_NED[1],
    a_NED[2] + g,
    q_dot.a,
    q_dot.b,
    q_dot.c,
    q_dot.d,
    #0,0,0,0,0,0
    wbx, wby, wbz, wbp, wbq, wbr
])

# discretized dynamics:
f = X + f_continuous*dt

# measurement model MOCAP:
use_quat = symbols('ekf_use_quat')
h_mocap = Matrix([x,y,z,use_quat*qw,use_quat*qx,use_quat*qy,use_quat*qz])

# output function BODY V (XY)
quat_inv = Quaternion(qw, -qx, -qy, -qz)
v_body = Quaternion.rotate_point([vx,vy,vz], quat_inv)
h_vbody = Matrix([v_body[0], v_body[1]])

# matrices:
F = f.jacobian(X)
L = f.jacobian(W)
H_mocap = h_mocap.jacobian(X)
H_vbody = h_vbody.jacobian(X)

# substitute W with 0
f = f.subs([(w,0) for w in W])
F = F.subs([(w,0) for w in W])
L = L.subs([(w,0) for w in W])
H_mocap = H_mocap.subs([(w,0) for w in W])
H_vbody = H_vbody.subs([(w,0) for w in W])

# extra matrices
symmetric_indexing = lambda i,j: int(j*(j+1)/2+i) if j>=i else int(i*(i+1)/2+j)
P = Matrix([[symbols(f'P[{symmetric_indexing(i,j)}]') for j in range(len(X))] for i in range(len(X))])
Q = Matrix([[symbols(f'Q[{i}]') if i==j else '0' for j in range(len(W))] for i in range(len(W))])
R_mocap = Matrix([[symbols(f'R_mocap[{i}]') if i==j else '0' for j in range(len(h_mocap))] for i in range(len(h_mocap))])
R_vbody = Matrix([[symbols(f'R_vbody[{i}]') if i==j else '0' for j in range(len(h_vbody))] for i in range(len(h_vbody))])
Z_mocap = Matrix([symbols(f'Z[{i}]') for i in range(len(h_mocap))])
Z_vbody = Matrix([symbols(f'Z[{i}]') for i in range(len(h_vbody))])

#%% generate code

# from sympy import symbols, sin
from sympy.codegen.ast import CodeBlock, Assignment

# EKF equations from: https://en.wikipedia.org/wiki/Extended_Kalman_filter: Non-additive noise formulation and equations

# PREDICTION STEP
Xpred = f
Ppred = F*P*F.T + L*Q*L.T

# S CALCULATION MOCAP
S_mocap = H_mocap*P*H_mocap.T + R_mocap
sdim_mocap = S_mocap.shape[0]
xdim = len(X)

# S CALCULATION BODYV
S_vbody = H_vbody*P*H_vbody.T + R_vbody
sdim_vbody = S_vbody.shape[0]

## new way with numerical inversion
HP_mocap = H_mocap*P
K_mocap = Matrix([[symbols(f'K_mocap[{j*xdim + i}]') for j in range(sdim_mocap)] for i in range(xdim)])

HP_vbody = H_vbody*P
K_vbody = Matrix([[symbols(f'K_vbody[{j*xdim + i}]') for j in range(sdim_vbody)] for i in range(xdim)])

Xup_mocap = X + K_mocap*(Z_mocap - h_mocap)
Pup_mocap = (eye(len(X)) - K_mocap*H_mocap)*P

Xup_vbody = X + K_vbody*(Z_vbody - h_vbody)
Pup_vbody = (eye(len(X)) - K_vbody*H_vbody)*P


# assignments
Xpred_assigments = [Assignment(symbols(f'X_new[{i}]'), Xpred[i]) for i in range(len(X))]
Ppred_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Ppred[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated

S_mocap_assigments = [Assignment(symbols(f'S_mocap[{j*sdim_mocap+i}]'), S_mocap[i,j]) for i in range(len(Z_mocap)) for j in range(len(Z_mocap))] # both triangular, full matrix
HP_mocap_assigments = [Assignment(symbols(f'HP_mocap[{j*sdim_mocap+i}]'), HP_mocap[i,j]) for j in range(xdim) for i in range(sdim_mocap)]
Xup_mocap_assigments = [Assignment(symbols(f'X_new[{i}]'), Xup_mocap[i]) for i in range(len(X))]
Pup_mocap_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Pup_mocap[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated

S_vbody_assigments = [Assignment(symbols(f'S_vbody[{j*sdim_vbody+i}]'), S_vbody[i,j]) for i in range(len(Z_vbody)) for j in range(len(Z_vbody))] # both triangular, full matrix
HP_vbody_assigments = [Assignment(symbols(f'HP_vbody[{j*sdim_vbody+i}]'), HP_vbody[i,j]) for j in range(xdim) for i in range(sdim_vbody)]
Xup_vbody_assigments = [Assignment(symbols(f'X_new[{i}]'), Xup_vbody[i]) for i in range(len(X))]
Pup_vbody_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Pup_vbody[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated


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

# UPDATE CODE FOR MOCAP-----------------------------------
# Prepare numerical solution for K
print('S and PHT MATRICES FOR SOLVING K:')
# code generation
s_code_mocap = CodeBlock(*S_mocap_assigments, *HP_mocap_assigments)
# common subexpression elimination with tmp variables
print('CSE')
s_code_mocap = s_code_mocap.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
# simplify
print('SIMPLIFY')
s_code_mocap = s_code_mocap.simplify()
# count number of tmp variables
s_code_mocap_N_tmps = len([s for s in s_code_mocap.left_hand_sides if s.name.startswith('tmp')])
# generate C code
print('CCODE')
s_code_mocap = ccode(s_code_mocap)


# UPDATE STEP
print('UPDATE:')
# code generation
update_code_mocap = CodeBlock(*Xup_mocap_assigments, *Pup_mocap_assigments)
# common subexpression elimination with tmp variables
print('CSE')
update_code_mocap = update_code_mocap.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
# simplify
print('SIMPLIFY')
update_code_mocap = update_code_mocap.simplify()
# count number of tmp variables
update_code_mocap_N_tmps = len([s for s in update_code_mocap.left_hand_sides if s.name.startswith('tmp')])
# generate C code
print('CCODE')
update_code_mocap = ccode(update_code_mocap)
# --------------------------------------------------------

# UPDATE CODE FOR VBODY-----------------------------------
# Prepare numerical solution for K
print('S and PHT MATRICES FOR SOLVING K:')
# code generation
s_code_vbody = CodeBlock(*S_vbody_assigments, *HP_vbody_assigments)
# common subexpression elimination with tmp variables
print('CSE')
s_code_vbody = s_code_vbody.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
# simplify
print('SIMPLIFY')
s_code_vbody = s_code_vbody.simplify()
# count number of tmp variables
s_code_vbody_N_tmps = len([s for s in s_code_vbody.left_hand_sides if s.name.startswith('tmp')])
# generate C code
print('CCODE')
s_code_vbody = ccode(s_code_vbody)


# UPDATE STEP
print('UPDATE:')
# code generation
update_code_vbody = CodeBlock(*Xup_vbody_assigments, *Pup_vbody_assigments)
# common subexpression elimination with tmp variables
print('CSE')
update_code_vbody = update_code_vbody.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
# simplify
print('SIMPLIFY')
update_code_vbody = update_code_vbody.simplify()
# count number of tmp variables
update_code_vbody_N_tmps = len([s for s in update_code_vbody.left_hand_sides if s.name.startswith('tmp')])
# generate C code
print('CCODE')
update_code_vbody = ccode(update_code_vbody)
# --------------------------------------------------------

#%% post-proc code

# replace sin, cos, tan, pow with sinf, cosf, tanf, powf
prediction_code = prediction_code.replace('sin(', 'sinf(')
prediction_code = prediction_code.replace('cos(', 'cosf(')
prediction_code = prediction_code.replace('tan(', 'tanf(')
prediction_code = prediction_code.replace('pow(', 'powf(')
prediction_code = prediction_code.replace('\n', '\n\t')
# mocap
s_code_mocap = s_code_mocap.replace('sin(', 'sinf(')
s_code_mocap = s_code_mocap.replace('cos(', 'cosf(')
s_code_mocap = s_code_mocap.replace('tan(', 'tanf(')
s_code_mocap = s_code_mocap.replace('pow(', 'powf(')
s_code_mocap = s_code_mocap.replace('\n', '\n\t')
update_code_mocap = update_code_mocap.replace('sin(', 'sinf(')
update_code_mocap = update_code_mocap.replace('cos(', 'cosf(')
update_code_mocap = update_code_mocap.replace('tan(', 'tanf(')
update_code_mocap = update_code_mocap.replace('pow(', 'powf(')
update_code_mocap = update_code_mocap.replace('\n', '\n\t')
# vbody
s_code_vbody = s_code_vbody.replace('sin(', 'sinf(')
s_code_vbody = s_code_vbody.replace('cos(', 'cosf(')
s_code_vbody = s_code_vbody.replace('tan(', 'tanf(')
s_code_vbody = s_code_vbody.replace('pow(', 'powf(')
s_code_vbody = s_code_vbody.replace('\n', '\n\t')
update_code_vbody = update_code_vbody.replace('sin(', 'sinf(')
update_code_vbody = update_code_vbody.replace('cos(', 'cosf(')
update_code_vbody = update_code_vbody.replace('tan(', 'tanf(')
update_code_vbody = update_code_vbody.replace('pow(', 'powf(')
update_code_vbody = update_code_vbody.replace('\n', '\n\t')


#%% emit code

from os import path, makedirs
from jinja2 import Environment, FileSystemLoader

# set data used in templates
data = {}
data['lenX'] = len(X)
data['lenQ'] = len(W)
data['lenU'] = len(U)
data['lenh_mocap'] = len(h_mocap)
data['lenh_vbody'] = len(h_vbody)
data['lenTmp'] = max(prediction_code_N_tmps, s_code_mocap_N_tmps, update_code_mocap_N_tmps, s_code_vbody_N_tmps, update_code_vbody_N_tmps)
data['prediction_code'] = prediction_code
data['prepare_gain_code_mocap'] = s_code_mocap
data['update_code_mocap'] = update_code_mocap
data['prepare_gain_code_vbody'] = s_code_vbody
data['update_code_vbody'] = update_code_vbody

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

