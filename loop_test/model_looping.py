#!/usr/bin/env python3

import sys
import argparse
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np


# need to figure out where to place loops, should be ok to place just before running solve_ivp

# Set copy numbers, degradation rates, and maximum protein production rates
# should there also be a production term for mRNA indicating promoter strength? or should copy number and promoter strength etc just be combined to 1 production rate for each gene?
production_m_OI, production_m_reg = 10, 2


# resource consumption per protein - should this be per protein produced?
res_consumption = 1


# contcentration constants -for now 
K_p_OI = 100000
K_p_reg = 1000
K_res = 500

# initial amounts of each product (mRNA and protein from each gene) and resource level (maybe number of ribosomes?)
m_OI_init = 10
p_OI_init = 1000
m_reg_init = 1
p_reg_init = 10
res_init = 100000

# Michaelis-menten-like (Hill functions, so far n = 1)
def m_m_positive(substrate_conc, K):
    return substrate_conc/(K+substrate_conc)

def m_m_negative(substrate_conc, K):
    return K/(K+substrate_conc)

# define stopping conditions
# if run out of resources or get negative values of protein
def run_out_res(t, y):
    return y[4]
run_out_res.terminal = True 

def run_out_p_OI(t, y):
    return y[1]
run_out_p_OI.terminal = True 

def run_out_p_reg(t, y):
    return y[3]
run_out_p_reg.terminal = True 

# write up equations
def model(indep, init_deps):
    '''
    Defining the set of differential equations to be solved
    Input: independent variable, initial values of dependent variables (list)
    Output: list of equations in the system 
    '''
    m_OI, p_OI, m_reg, p_reg, resource = init_deps

    # equations
    dm_OI = m_m_negative(p_reg, K_p_reg) * copy_num_1 - deg_rate_m * m_OI
    dp_OI = max_prod * m_m_positive(resource, K_res) * m_OI - deg_rate_p * p_OI
    dm_reg = m_m_positive(p_OI, K_p_OI) * copy_num_2 - deg_rate_m * m_reg
    dp_reg = max_prod * m_m_positive(resource, K_res) * m_reg - deg_rate_p * p_reg
    # assume constant resource consumption for protein for now
    dres = - res_consump * (p_OI + p_reg) # + something to replenish ? + res_init?

    return [dm_OI, dp_OI, dm_reg, dp_reg, dres]

# list of initial conditions
inits = [m_OI_init, p_OI_init, m_reg_init, p_reg_init, res_init]

# time steps
t = np.linspace(0, 20, 100)

# parameter spaces
deg_rates_m = np.linspace(0.1, 100, 20)
deg_rates_p = np.linspace(1, 100, 20)
max_prods = np.linspace(100, 10000, 20)

for deg_rate_m in range(11):
    for deg_rate_p in range(10,101,10):
        for max_prod in range(10,10001,1000):        
            # solving the ODE system
            solution = ig.solve_ivp(model, [min(t),max(t)], inits, t_eval = t, events = [run_out_res, run_out_p_OI,run_out_p_reg])
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
            plt.title('resource level')
            plt.savefig(fname=f'mRNA_deg_{deg_rate_m}_prot_deg_{deg_rate_p}_maximum prod_{max_prod}.png')
            plt.clf()
#plt.legend(['m_OI', 'm_reg', 'p_reg', 'fitness index'])
