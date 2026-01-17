from sympy import *
from sympy.codegen.ast import CodeBlock, Assignment

# states, inputs, and process noise
X=Matrix(symbols('X[0] X[1] X[2] X[3] X[4] X[5] X[6] X[7] X[8] X[9] X[10] X[11] X[12] X[13] X[14] X[15] X[16] X[17] X[18] X[19] X[20] X[21] X[22] X[23] X[24] X[25]'))
U=Matrix(symbols('U[0] U[1] U[2] U[3] U[4] U[5]'))
W=Matrix(symbols('W[0] W[1] W[2] W[3] W[4] W[5] W[6] W[7] W[8] W[9] W[10] W[11] W[12] W[13] W[14] W[15] W[16] W[17] W[18] W[19] W[20] W[21]'))

# time step
dt=symbols('dt')

# continuous dynamics:
g = 9.81
# pos, vel, att, acc bias, gyro bias, extrinsics
x,y,z,vbx,vby,vbz,qw,qx,qy,qz,bx,by,bz,bp,bq,br,ex,ey,ez,ephi,etheta,epsi,fx,fy,cx,cy = X
ax,ay,az,p,q,r = U
wx,wy,wz,wp,wq,wr,wbx,wby,wbz,wbp,wbq,wbr,wex,wey,wez,wephi,wetheta,wepsi,wfx,wfy,wcx,wcy = W

quat = Quaternion(qw, qx, qy, qz) # does norm 1 automatically renormaize?
quat_inv = Quaternion(qw, -qx, -qy, -qz)
Rx = Matrix([[1, 0, 0], [0, cos(ephi), -sin(ephi)], [0, sin(ephi), cos(ephi)]])
Ry = Matrix([[cos(etheta), 0, sin(etheta)],[0, 1, 0],[-sin(etheta), 0, cos(etheta)]])
Rz = Matrix([[cos(epsi), -sin(epsi), 0],[sin(epsi), cos(epsi), 0], [0, 0, 1]])
Re = Rz*Ry*Rx

acc = Matrix([ax-bx-wx, ay-by-wy, az-bz-wz])
omega = Matrix([p-bp-wp, q-bq-wq, r-br-wr])

gbody = Quaternion.rotate_point([0, 0, g], quat_inv)
v_world = Quaternion.rotate_point([vbx, vby, vbz], quat)
v_body = Matrix([vbx, vby, vbz])

pqr_hat = Quaternion(0, p-bp-wp, q-bq-wq, r-br-wr)
q_dot = 0.5 * quat * pqr_hat

v_body_dot = acc - omega.cross(v_body) + Matrix(gbody)

f_continuous = Matrix([
    v_world[0], v_world[1], v_world[2],
    v_body_dot[0], v_body_dot[1], v_body_dot[2],
    q_dot.a,
    q_dot.b,
    q_dot.c,
    q_dot.d,
    wbx, wby, wbz, wbp, wbq, wbr,
    wex, wey, wez, wephi, wetheta, wepsi,
    wfx, wfy, wcx, wcy
])

# discretized dynamics:
f = X + f_continuous*dt

# NORMALIZATION:
q_norm = sqrt(qw**2 + qx**2 + qy**2 + qz**2)
# clamp
def clamp(v, vmin, vmax): return Max(Min(v, vmax), vmin)

# bias bounds
bx_min, bx_max, by_min, by_max, bz_min, bz_max, bp_min, bp_max, bq_min, bq_max, br_min, br_max = symbols('bx_min bx_max by_min by_max bz_min bz_max bp_min bp_max bq_min bq_max br_min br_max')
# extrinsic bounds
ex_min, ex_max, ey_min, ey_max, ez_min, ez_max, ephi_min, ephi_max, etheta_min, etheta_max, epsi_min, epsi_max = symbols('ex_min ex_max ey_min ey_max ez_min ez_max ephi_min ephi_max etheta_min etheta_max epsi_min epsi_max')
# intrinsics bounds
fx_min, fx_max, fy_min, fy_max, cx_min, cx_max, cy_min, cy_max = symbols('fx_min fx_max fy_min fy_max cx_min cx_max cy_min cy_max')

Xnorm = Matrix([
    x, y, z,
    vbx, vby, vbz,
    qw/q_norm, qx/q_norm, qy/q_norm, qz/q_norm,
    clamp(bx, bx_min, bx_max),
    clamp(by, by_min, by_max),
    clamp(bz, bz_min, bz_max),
    clamp(bp, bp_min, bp_max),
    clamp(bq, bq_min, bq_max),
    clamp(br, br_min, br_max),
    clamp(ex, ex_min, ex_max),
    clamp(ey, ey_min, ey_max),
    clamp(ez, ez_min, ez_max),
    clamp(ephi, ephi_min, ephi_max),
    clamp(etheta, etheta_min, etheta_max),
    clamp(epsi, epsi_min, epsi_max),
    clamp(fx, fx_min, fx_max),
    clamp(fy, fy_min, fy_max),
    clamp(cx, cx_min, cx_max),
    clamp(cy, cy_min, cy_max)
])

# collect all min/max bounds in a single global_vars list
global_vars = [
    bx_min, bx_max, by_min, by_max, bz_min, bz_max,
    bp_min, bp_max, bq_min, bq_max, br_min, br_max,
    ex_min, ex_max, ey_min, ey_max, ez_min, ez_max,
    ephi_min, ephi_max, etheta_min, etheta_max, epsi_min, epsi_max,
    fx_min, fx_max, fy_min, fy_max, cx_min, cx_max, cy_min, cy_max
]
global_vars_default = [-10,10]*12
# intrinsics reasonable defaults
global_vars_default += [0,2000,0,2000,0,2000, 0,2000]
global_vars = list(zip([v.name for v in global_vars], global_vars_default))

# MEASUREMENT MODELS
measurement_models = {}

# position measurement
h = Matrix([x,y,z])
H = h.jacobian(X)
measurement_models['pos'] = (h, H)

# position + quat measurement
h = Matrix([x,y,z,qw,qx,qy,qz])
H = h.jacobian(X)
measurement_models['pos_quat'] = (h, H)

# body velocity measurement
h = Matrix([v_body[0], v_body[1], v_body[2]])
H = h.jacobian(X)
measurement_models['vel_body'] = (h, H)

# projected points measurement 1-16
max_points = 16
points_3d = [symbols(f'p3d{i}_x p3d{i}_y p3d{i}_z') for i in range(max_points)]
points_2d = []
for p3d in points_3d:
    # point
    px, py, pz = p3d
    # from world to body coordinates
    p_body = Quaternion.rotate_point([px-x, py-y, pz-z], quat_inv)
    px, py, pz = p_body
    # from body to camera coordinates
    # p_cam = Quaternion.rotate_point([px-ex, py-ey, pz-ez], Quaternion(eqw, -eqx, -eqy, -eqz))
    p_cam = Re.T * Matrix([px-ex, py-ey, pz-ez])
    px, py, pz = p_cam
    # camera coordinates to opencv convention
    px, py, pz = py, pz, px
    # project to 2d
    u = px*fx/pz + cx
    v = py*fy/pz + cy
    points_2d.append((u,v))

h_ = []
for p in points_2d:
    h_.append(p[0])
    h_.append(p[1])
# h = Matrix(h)

for i in range(max_points):
    hi = Matrix(h_[0:2*(i+1)])
    Hi = hi.jacobian(X)
    params = []
    for p in points_3d[0:i+1]:
        params.extend(p)
    measurement_models[f'points_{i+1}'] = (hi, Hi, params)

# CODE GENERATION
def renaming(code):
    code = code.replace('sin(', 'sinf(')
    code = code.replace('cos(', 'cosf(')
    code = code.replace('tan(', 'tanf(')
    code = code.replace('pow(', 'powf(')
    # code = code.replace('\n', '\n\t')
    return code

symmetric_indexing = lambda i,j: int(j*(j+1)/2+i) if j>=i else int(i*(i+1)/2+j)
P = Matrix([[symbols(f'P[{symmetric_indexing(i,j)}]') for j in range(len(X))] for i in range(len(X))])
Q = Matrix([[symbols(f'Q[{i}]') if i==j else '0' for j in range(len(W))] for i in range(len(W))])

# Jacobians
F = f.jacobian(X)
L = f.jacobian(W)

# substitute W with 0
f = f.subs([(w,0) for w in W])
F = F.subs([(w,0) for w in W])
L = L.subs([(w,0) for w in W])

Xpred = f
Ppred = F*P*F.T + L*Q*L.T

# tmp variables
num_tmps = 0

# assignments
Xpred_assigments = [Assignment(symbols(f'X_new[{i}]'), Xpred[i]) for i in range(len(X))]
Ppred_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Ppred[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated

# code generation
prediction_code = CodeBlock(*Xpred_assigments, *Ppred_assigments)
prediction_code = prediction_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
prediction_code = prediction_code.simplify()
prediction_code_N_tmps = len([s for s in prediction_code.left_hand_sides if s.name.startswith('tmp')])
num_tmps = max(num_tmps, prediction_code_N_tmps)
prediction_code = ccode(prediction_code)
prediction_code = renaming(prediction_code)

# integrate code
# X integration code is just like prediction code the X is renamed to X_
Xint = Xpred.subs([(symbols(f'X[{i}]'), symbols(f'X_[{i}]')) for i in range(len(X))])
Xint_assignments = [Assignment(symbols(f'X_integrate[{i}]'), Xint[i]) for i in range(len(X))]
integrate_code = CodeBlock(*Xint_assignments)
integrate_code = integrate_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
integrate_code = integrate_code.simplify()
integrate_code_N_tmps = len([s for s in integrate_code.left_hand_sides if s.name.startswith('tmp')])
num_tmps = max(num_tmps, integrate_code_N_tmps)
integrate_code = ccode(integrate_code)
integrate_code = renaming(integrate_code)

# normalization code
Xnorm_assigments = [Assignment(symbols(f'X_new[{i}]'), Xnorm[i]) for i in range(len(X))]
norm_code = CodeBlock(*Xnorm_assigments)
norm_code = norm_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
norm_code = norm_code.simplify()
norm_code_N_tmps = len([s for s in norm_code.left_hand_sides if s.name.startswith('tmp')])
num_tmps = max(num_tmps, norm_code_N_tmps)
norm_code = ccode(norm_code)
norm_code = renaming(norm_code)

# UPDATE CODE
measurements = []
for name, val in measurement_models.items():
    h, H = val[0], val[1]
    if len(val) > 2:
        params = val[2]
    else:
        params = []
    print(params)
    
    num_measurements = len(h)
    num_states = len(X)
    print(f'Generating code for measurement model {name} with {num_measurements} measurements')
    
    # symbols
    Z = Matrix([symbols(f'Z[{i}]') for i in range(num_measurements)])
    R = Matrix([[symbols(f'R[{i}]') if i==j else '0' for j in range(num_measurements)] for i in range(num_measurements)])
    K = Matrix([[symbols(f'K[{j*num_states + i}]') for j in range(num_measurements)] for i in range(num_states)])
    
    # S = H*P*H.T + R
    # HP = H*P
    # sdim = S.shape[0]
    # xdim = len(X)
    
    # Xup = X + K*(Z - h)
    # Pup = (eye(len(X)) - K*H)*P
    # joseph form to ensure positive semi-definiteness
    # Pup = (eye(len(X)) - K*H)*P*(eye(len(X)) - K*H).T + K*R*K.T
    
    # S_assigments = [Assignment(symbols(f'S[{j*sdim+i}]'), S[i,j]) for i in range(len(Z)) for j in range(len(Z))] # both triangular, full matrix
    # HP_assigments = [Assignment(symbols(f'HP[{j*sdim+i}]'), HP[i,j]) for j in range(xdim) for i in range(sdim)]

    # Xup_assigments = [Assignment(symbols(f'X_new[{i}]'), Xup[i]) for i in range(len(X))]
    # Pup_assigments = [Assignment(symbols(f'P_new[{symmetric_indexing(i,j)}]'), Pup[i,j]) for i in range(len(X)) for j in range(len(X)) if j >= i] # only the lower diagonal will be calculated

    # code generation
    # s_code = CodeBlock(*S_assigments, *HP_assigments)
    # s_code = s_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
    # s_code = s_code.simplify()
    # s_code_N_tmps = len([s for s in s_code.left_hand_sides if s.name.startswith('tmp')])
    # num_tmps = max(num_tmps, s_code_N_tmps)
    # s_code = ccode(s_code)
    # s_code = renaming(s_code)
    
    # update_X_code = CodeBlock(*Xup_assigments)
    # update_X_code = update_X_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
    # update_X_code = update_X_code.simplify()
    # update_X_code_N_tmps = len([s for s in update_X_code.left_hand_sides if s.name.startswith('tmp')])
    # num_tmps = max(num_tmps, update_X_code_N_tmps)
    # update_X_code = ccode(update_X_code)
    # update_X_code = renaming(update_X_code)
    
    # update_P_code = CodeBlock(*Pup_assigments)
    # update_P_code = update_P_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
    # update_P_code = update_P_code.simplify()
    # update_P_code_N_tmps = len([s for s in update_P_code.left_hand_sides if s.name.startswith('tmp')])
    # num_tmps = max(num_tmps, update_P_code_N_tmps)
    # update_P_code = ccode(update_P_code)
    # update_P_code = renaming(update_P_code)
    
    h_assignments = [Assignment(symbols(f'h[{i}]'), h[i]) for i in range(len(h))]
    h_code = CodeBlock(*h_assignments)
    h_code = h_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
    h_code = h_code.simplify()
    h_code_N_tmps = len([s for s in h_code.left_hand_sides if s.name.startswith('tmp')])
    num_tmps = max(num_tmps, h_code_N_tmps)
    h_code = ccode(h_code)
    h_code = renaming(h_code)
    
    H_assignments = [Assignment(symbols(f'H[{i + j*len(h)}]'), H[i,j]) for i in range(len(h)) for j in range(len(X))]
    H_code = CodeBlock(*H_assignments)
    H_code = H_code.cse(symbols=(symbols(f'tmp[{i}]') for i in range(10000)))
    H_code = H_code.simplify()
    H_code_N_tmps = len([s for s in H_code.left_hand_sides if s.name.startswith('tmp')])
    num_tmps = max(num_tmps, H_code_N_tmps)
    H_code = ccode(H_code)
    H_code = renaming(H_code)
    
    # update_code[name] = (s_code, u_code, max(s_code_N_tmps, u_code_N_tmps))
    measurements.append({
        'name': name,
        'num_measurements': num_measurements,
        'num_params': len(params),
        'param_names': [p.name for p in params],
        'h_code': h_code,
        'H_code': H_code
    })
    
# fill in the jinja template    
from os import path, makedirs
from jinja2 import Environment, FileSystemLoader

# set data used in templates
data = {}
data['lenX'] = len(X)
data['lenQ'] = len(W)
data['lenU'] = len(U)
data['prediction_code'] = prediction_code
data['integrate_code'] = integrate_code
data['norm_code'] = norm_code
data['global_vars'] = global_vars
data['measurements'] = measurements
data['num_tmps'] = num_tmps

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