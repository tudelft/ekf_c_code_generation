from sympy import *

# states, inputs, and process noise
X=Matrix(symbols('X[0] X[1] X[2] X[3] X[4] X[5] X[6] X[7] X[8] X[9] X[10] X[11] X[12] X[13] X[14] X[15]'))
U=Matrix(symbols('U[0] U[1] U[2] U[3] U[4] U[5]'))
W=Matrix(symbols('W[0] W[1] W[2] W[3] W[4] W[5] W[6] W[7] W[8] W[9] W[10] W[11]'))

# time step
dt=symbols('dt')

g = 9.81
# state
x,y,z,vx,vy,vz,qw,qx,qy,qz,lx,ly,lz,lp,lq,lr = X
# input
ax,ay,az,p,q,r = U
# process noise
wx,wy,wz,wp,wq,wr,wbx,wby,wbz,wbp,wbq,wbr = W

quat = Quaternion(qw, qx, qy, qz) # does norm 1 automatically renormaize?
a_NED = Quaternion.rotate_point([ax-lx-wx,ay-ly-wy,az-lz-wz], quat)
quat_inv = Quaternion(qw, -qx, -qy, -qz)
v_body = Quaternion.rotate_point([vx,vy,vz], quat_inv)
vbx, vby, vbz = v_body

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

phi = atan2(2*(qw*qx + qy*qz), 1 - 2*(qx**2 + qy**2))
theta = asin(2*(qw*qy - qz*qx))
psi = atan2(2*(qw*qz + qx*qy), 1 - 2*(qy**2 + qz**2))

# state transition model
f = Matrix([
    x+vx*dt,
    y+vy*dt,
    z+vz*dt,
    vx + a_NED[0]*dt,
    vy + a_NED[1]*dt,
    vz + (a_NED[2]+g)*dt,
    qw + q_dot.a*dt,
    qx + q_dot.b*dt,
    qy + q_dot.c*dt,
    qz + q_dot.d*dt,
    lx + wbx,   # acc x bias
    ly + wby,   # acc y bias
    lz + wbz,   # acc z bias
    lp + wbp,   # gyro x bias
    lq + wbq,   # gyro y bias
    lr + wbr,   # gyro z bias
])

# output function (measurement model):
h_pnp = Matrix([x,y,z,qw,qx,qy,qz])
h_v_body = Matrix([vx,vy,vz,phi,theta])

# matrices:
F = f.jacobian(X)
L = f.jacobian(W)
H_pnp = h_pnp.jacobian(X)
H_v_body = h_v_body.jacobian(X)

# substitute W with 0
f = f.subs([(w,0) for w in W])
F = F.subs([(w,0) for w in W])
L = L.subs([(w,0) for w in W])
H_pnp = H_pnp.subs([(w,0) for w in W])
H_v_body = H_v_body.subs([(w,0) for w in W])

# extra matrices
symmetric_indexing = lambda i,j: int(j*(j+1)/2+i) if j>=i else int(i*(i+1)/2+j)
P = Matrix([[symbols(f'P[{symmetric_indexing(i,j)}]') for j in range(len(X))] for i in range(len(X))])
Q = Matrix([[symbols(f'Q[{i}]') if i==j else '0' for j in range(len(W))] for i in range(len(W))])
R_pnp = Matrix([[symbols(f'R_pnp[{i}]') if i==j else '0' for j in range(len(h_pnp))] for i in range(len(h_pnp))])
R_v_body = Matrix([[symbols(f'R_v_body[{i}]') if i==j else '0' for j in range(len(h_v_body))] for i in range(len(h_v_body))])

Z_pnp = Matrix([symbols(f'Z_pnp[{i}]') for i in range(len(h_pnp))])
Z_v_body = Matrix([symbols(f'Z_v_body[{i}]') for i in range(len(h_v_body))])


#%% generate code

# from sympy import symbols, sin
from sympy.codegen.ast import CodeBlock, Assignment

# EKF equations from: https://en.wikipedia.org/wiki/Extended_Kalman_filter: Non-additive noise formulation and equations

# PREDICTION STEP
Xpred = f
Ppred = F*P*F.T + L*Q*L.T

# assignments
Xpred_assigments = [Assignment(symbols(f'X_new[{i}]'), Xpred[i]) for i in range(len(X))]
Ppred_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Ppred[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated

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

def get_update_code(h, H, R, Z, name=''):
    # UPDATE STEP
    S = H*P*H.T + R
    sdim = S.shape[0]
    xdim = len(X)
    HP = H*P
    K = Matrix([[symbols(f'K_{name}[{j*xdim + i}]') for j in range(sdim)] for i in range(xdim)])
    Xup = X + K*(Z - h)
    Pup = (eye(len(X)) - K*H)*P
    
    S_assigments = [Assignment(symbols(f'S_{name}[{j*sdim+i}]'), S[i,j]) for i in range(len(Z)) for j in range(len(Z))] # both triangular, full matrix
    HP_assigments = [Assignment(symbols(f'HP_{name}[{j*sdim+i}]'), HP[i,j]) for j in range(xdim) for i in range(sdim)]
    Xup_assigments = [Assignment(symbols(f'X_new[{i}]'), Xup[i]) for i in range(len(X))]
    Pup_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Pup[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated

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
    
    return update_code, s_code, update_code_N_tmps, s_code_N_tmps
    

# PNP UPDATE CODE
update_code_pnp, s_code_pnp, update_code_N_tmps_pnp, s_code_N_tmps_pnp = get_update_code(h_pnp, H_pnp, R_pnp, Z_pnp, name='pnp')
# V_BODY UPDATE CODE
update_code_v_body, s_code_v_body, update_code_N_tmps_v_body, s_code_N_tmps_v_body = get_update_code(h_v_body, H_v_body, R_v_body, Z_v_body, name='v_body')


#%% post-proc code

# replace sin, cos, tan, pow with sinf, cosf, tanf, powf
def float_functions(txt):
    txt = txt.replace('sin(', 'sinf(')
    txt = txt.replace('cos(', 'cosf(')
    txt = txt.replace('tan(', 'tanf(')
    txt = txt.replace('pow(', 'powf(')
    txt = txt.replace('\n', '\n\t')
    return txt

# prediction code
prediction_code = float_functions(prediction_code)
# pnp update code
s_code_pnp = float_functions(s_code_pnp)
update_code_pnp = float_functions(update_code_pnp)
# acc update code
s_code_v_body = float_functions(s_code_v_body)
update_code_v_body = float_functions(update_code_v_body)

#%% emit code

from os import path, makedirs
from jinja2 import Environment, FileSystemLoader

# set data used in templates
data = {}
data['lenX'] = len(X)
data['lenQ'] = len(W)
data['lenU'] = len(U)
data['lenh_pnp'] = len(h_pnp)
data['lenh_v_body'] = len(h_v_body)
data['lenTmp'] = max(prediction_code_N_tmps, s_code_N_tmps_pnp, update_code_N_tmps_pnp, s_code_N_tmps_v_body, update_code_N_tmps_v_body)
data['prediction_code'] = prediction_code
data['prepare_gain_code_pnp'] = s_code_pnp
data['update_code_pnp'] = update_code_pnp
data['prepare_gain_code_v_body'] = s_code_v_body
data['update_code_v_body'] = update_code_v_body

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

