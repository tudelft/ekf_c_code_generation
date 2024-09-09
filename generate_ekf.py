from sympy import *

# states, inputs, and process noise
X=Matrix(symbols('X[0] X[1] X[2] X[3] X[4] X[5] X[6] X[7] X[8] X[9] X[10] X[11] X[12] X[13] X[14]'))
U=Matrix(symbols('U[0] U[1] U[2] U[3] U[4] U[5]'))
W=Matrix(symbols('W[0] W[1] W[2] W[3] W[4] W[5]'))

# time step
dt=symbols('dt')

# continuous dynamics:
g = 9.81
x,y,z,vx,vy,vz,phi,theta,psi,lx,ly,lz,lp,lq,lr = X
ax,ay,az,p,q,r = U
wx,wy,wz,wp,wq,wr = W

Rx = Matrix([[1, 0, 0], [0, cos(phi), -sin(phi)], [0, sin(phi), cos(phi)]])
Ry = Matrix([[cos(theta), 0, sin(theta)],[0, 1, 0],[-sin(theta), 0, cos(theta)]])
Rz = Matrix([[cos(psi), -sin(psi), 0],[sin(psi), cos(psi), 0], [0, 0, 1]])
R = Rz*Ry*Rx
a_NED = R*Matrix([ax-lx-wx,ay-ly-wy,az-lz-wz])

f_continuous = Matrix([
    vx,
    vy,
    vz,
    a_NED[0],
    a_NED[1],
    a_NED[2] + g,
    (p-lp-wp) + (q-lq-wq)*sin(phi)*tan(theta) + (r-lr-wr)*cos(phi)*tan(theta),
    (q-lq-wq)*cos(phi) - (r-lr-wr)*sin(phi),
    (q-lq-wq)*sin(phi)/cos(theta) + (r-lr-wr)*cos(phi)/cos(theta),
    0,0,0,0,0,0
])


# discretized dynamics:
f = X + f_continuous*dt

# output function (measurement model):
use_phi, use_theta, use_psi = symbols('ekf_use_phi ekf_use_theta ekf_use_psi') # 0 or 1
h = Matrix([x,y,z,use_phi*phi,use_theta*theta,use_psi*psi])

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



#%% generate code

# from sympy import symbols, sin
from sympy.codegen.ast import CodeBlock, Assignment

# function for inverting matrix with symbolic elements that is quicker
def invert_function(num):
    M = Matrix([[symbols(f'M[{i}{j}]') for j in range(num)] for i in range(num)])
    M_inv = M.inverse_LU()
    return lambda M_: M_inv.subs({M[i,j]:M_[i,j] for i in range(num) for j in range(num)})


# EKF equations from: https://en.wikipedia.org/wiki/Extended_Kalman_filter: Non-additive noise formulation and equations

# PREDICTION STEP
Xpred = f
Ppred = F*P*F.T + L*Q*L.T

# UPDATE STEP
S = H*P*H.T + R
inv__ = invert_function(S.shape[0])
K = P*H.T*inv__(S)
Xup = X + K*(Z - h)
Pup = (eye(len(X)) - K*H)*P


# assignments
Xpred_assigments = [Assignment(symbols(f'X_new[{i}]'), Xpred[i]) for i in range(len(X))]
Ppred_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Ppred[i,j]) for i in range(len(X)) for j in range(len(X)) if i >= j] # only the lower diagonal will be calculated

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
update_code = update_code.replace('sin(', 'sinf(')
update_code = update_code.replace('cos(', 'cosf(')
update_code = update_code.replace('tan(', 'tanf(')
update_code = update_code.replace('pow(', 'powf(')
update_code = update_code.replace('\n', '\n\t')

#print('Prediction code:')
#print(prediction_code)
#print('Update code:')
#print(update_code)


#%% emit code

from os import path, makedirs
from jinja2 import Environment, FileSystemLoader

# set data used in templates
data = {}
data['lenX'] = len(X)
data['lenU'] = len(U)
data['lenh'] = len(h)
data['lenTmp'] = max(prediction_code_N_tmps, update_code_N_tmps)
data['prediction_code'] = prediction_code
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

