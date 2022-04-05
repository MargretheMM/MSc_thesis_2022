#!/usr/bin/env python3

import sys
import argparse
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np

# Set degradation rates and maximum production rates for mRNA and protein
deg_rate_m = 10
deg_rate_p = 1
production_m_v, production_m_r = 100, 5
max_protein_production = 1000



''' Resource terms may all need to be re-framed into fractions of maximum possible resources?'''
# resource consumption per protein - should this be per protein produced?
store_consumption = 1
# max resource level
max_store = 100000
# resource replenishment rate
store_replenish = 1000
# resource dilution rate - due to growth
store_dilution = 10



# contcentration constants -for now 
K_p_v = 100000
K_p_r = 1000
K_store = 500

# initial amounts of each product (mRNA and protein from each gene) and resource level (maybe number of ribosomes?)
m_v_init = 10
p_v_init = 100
m_r_init = 1
p_r_init = 10
store_init = 80000

# Control functions. Currently Michaelis-menten-like (Hill functions, so far n = 1)
def control_positive(substrate_conc, K):
    return substrate_conc/(K+substrate_conc)

def control_negative(substrate_conc, K):
    return K/(K+substrate_conc)

# define stopping conditions
# if run out of resources or get negative values of protein
def run_out_store(t, y):
    return y[4]
run_out_store.terminal = True 

def run_out_p_v(t, y):
    return y[1]
run_out_p_v.terminal = True 

def run_out_p_r(t, y):
    return y[3]
run_out_p_r.terminal = True 

# write up equations
def model(indep, init_deps):
    '''
    Defining the set of differential equations to be solved
    Input: independent variable, initial values of dependent variables (list)
    Output: list of equations in the system 
    '''
    m_v, p_v, m_r, p_r, store = init_deps

    # equations
    dm_v = control_negative(p_r, K_p_r) * production_m_v - deg_rate_m * m_v
    dp_v = max_protein_production * control_positive(store, K_store) * m_v - deg_rate_p * p_v
    dm_r = control_positive(p_v, K_p_v) * production_m_r - deg_rate_m * m_r
    dp_r = max_protein_production * control_positive(store, K_store) * m_r - deg_rate_p * p_r
    # assume constant resource consumption for protein for now
    dstore = store_replenish - store_consumption * p_v - store_dilution

    return [dm_v, dp_v, dm_r, dp_r, dstore]

# list of initial conditions
inits = [m_v_init, p_v_init, m_r_init, p_r_init, store_init]

# time steps
t = np.linspace(0, 20, 100)

# solving the ODE system
solution = ig.solve_ivp(model, [min(t),max(t)], inits, t_eval = t, events = [run_out_store, run_out_p_v,run_out_p_r])
plt.figure(1)
plt.subplot(511)
plt.plot(solution.t, solution.y[1].T)
plt.title('POI over time')
plt.subplot(512)
plt.plot(solution.t, solution.y[0].T)
plt.title('mRNA of gene 1')
plt.subplot(513)
plt.plot(solution.t, solution.y[2].T)
plt.title('mRNA of gene 2')
plt.subplot(514)
plt.plot(solution.t, solution.y[3].T)
plt.title('regulator protein')
plt.subplot(515)
plt.plot(solution.t, solution.y[4].T)
plt.title('store level')
#plt.legend(['m1', 'm2', 'p2', 'fitness index'])
plt.show()
