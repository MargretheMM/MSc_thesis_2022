#!/usr/bin/env python3

import numpy as np
import sympy as sp

# Model parameters - biologically based
R_productive = 1.5e5
p_max = 3e7
m_max = 3e4
beta_m = 1/4
beta_p = 0.05
beta_R = 2000/60
alpha_m = np.log(2)/1200
alpha_p = np.log(2)/18000
alpha_R = np.log(2)/21600
k = np.log(2)/6000
K_m_v = 1e3
K_R = 4.5e4
m_v_init = 0.05
p_v_init = 0.1
m_r_init = 0
p_r_init = 0
R_init = 1.2e5 

# Functions and model
# Control functions. Currently Michaelis-menten-like (Hill functions, so far n = 1)
def control_positive(K, substrate):
    return substrate/(K+substrate)

def control_negative(K, substrate):
    return K/(K+substrate)


# defining symbols
m_v, p_v, m_r, p_r, R = sp.symbols('m_v, p_v, m_r, p_r, R', real = True)
delta_m_v, delta_m_r, multi_m_r, multi_p_r, K_p_v, K_p_r, K_m_v, K_R = sp.symbols('delta_m_v, delta_m_r, multi_m_r, multi_p_r, K_p_v, K_p_r, K_m_v, K_R')
beta_p, beta_m, beta_R, alpha_R, alpha_m_v, alpha_p_v, m_max, R_productive, k = sp.symbols('beta_p, beta_m, beta_R, alpha_R, alpha_m_v, alpha_p_v, m_max, R_productive, k')


# equilibrium search equations

eq1 = sp.Eq(beta_m * delta_m_v *  control_negative(2, 2)  - alpha_m_v * m_max,0)
#print(sp.solve(eq1,delta_m_v))

eq2 = sp.Eq(beta_p * R * m_v/m_max * m_v * control_positive(K_R, R) * control_positive(K_m_v,m_v) - (alpha_p_v + k * (1 - m_v/m_max) ) * p_v,0)
#print(sp.solve(eq2,m_v))

eq3 = sp.Eq(beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_v * multi_m_r * m_r,0)
#print(sp.solve(eq3,[delta_m_r, multi_m_r]))

eq4 = sp.Eq(beta_p * 0.5 * R_productive * (1 - 1e3/m_max) * 0.5/(m_max - 1e3)  - (alpha_p_v * multi_p_r + k * (1 - 1e3/m_max)) * 2,0)
#print(sp.solve(eq4,multi_p_r))


eq5 = sp.Eq(beta_R * R * (1 - m_v/m_max) / R_productive - (alpha_R + k * R * (1 - m_v/m_max)/(R_productive)) * R,0)
print(sp.solve(eq5,R))

# Jacobian of model equation matrix
# equation matrix
eqns = sp.Matrix([beta_m * delta_m_v *  control_negative(K_p_r, p_r)  - alpha_m_v * m_v, beta_p * R * m_v/m_max * m_v * control_positive(K_R,R) * control_positive(K_m_v,m_v) - (alpha_p_v + k * (1 - m_v/m_max) ) * p_v, beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_v * multi_m_r * m_r, beta_p * m_r * R * (1 - m_v/m_max) * m_r/(m_max - m_v)  - (alpha_p_v * multi_p_r + k * (1 - m_v/m_max)) * p_r, beta_R * R * (1 - m_v/m_max) / R_productive - (alpha_R + k * (1 - m_v/m_max)) * R])
#variables
variables = sp.Matrix([m_v, p_v, m_r, p_r, R])
#jacobian
eqns.jacobian(variables)
