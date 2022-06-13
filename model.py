#!/usr/bin/env python3

import sys
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np


''' Model parameters - biologically based '''

# maximum number of ribosomes available for protein production - currently assume 95% of total - Metzl-Raz 2017 and von der Haar 2008
R= 2e5*0.95


# total number of RNA polymerases in cell - Borggrefe 2001
total_pol = 30000

# Maximum number of protein bound amino acids in cell - futcher et al 1999
p_max = 3e10

# Maximum transcription and translation rates
# max transcription in codons per hour per cell : codons/sec * sec/hour *number of polymerases = codons/hour - Yu and Nielsen 2019
beta_m = 50/3 * 3600 * total_pol
# maximum translation in aas per hour per cell : aa/ribosome/sec * sec/hour * number of ribsomes which can be diverted = aa/hour
delta_R = 10 * 3600
beta_p = delta_R * R


# Set degradation rates and maximum production rates for mRNA and protein - mRNA: Curran et al 2013, Perez-Ortin et al 2007, protein Christiano 2020
# rates based on half-life in h
alpha_m_v = np.log(2)/0.33
alpha_p_v = np.log(2)/5


''' Presumed engineerable paramters '''
# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 1e6
K_p_r = 50

# mRNA production rates - mix of copy number, promotor strength etc - will be scaled with transcription rate
delta_m_v, delta_m_r = 1, 1e-4

# regulator degradation rates
alpha_m_r = np.log(2)/0.33
alpha_p_r = np.log(2)/1

''' Initial conditions'''
# initial amounts of each product (mRNA and protein from each gene)
m_v_init = 5
p_v_init = 100
m_r_init = 1
p_r_init = 0
R_v_init = 5
R_o_init = R-R_v_init


'''Functions and model '''
# Control functions. Currently Michaelis-menten-like (Hill functions, so far n = 1)
def control_positive(K, substrate):
    return substrate/(K+substrate)

def control_negative(K, substrate):
    return K/(K+substrate)


# write up equations
def model(indep: float, init_deps):
    '''
    Defining the set of differential equations to be solved
    Input: independent variable, initial values of dependent variables (list)
    Output: list of equations in the system 
    '''
    m_v, p_v, m_r, p_r, R_v, R_o = init_deps
    
    # protein level in model
    p_tot = p_v + p_r + R * 12e3
    
    # equations
    dm_v = beta_m * delta_m_v *  control_negative(K_p_r, p_r)  - alpha_m_v * m_v
    
    if (0 < p_tot < p_max):
        dp_v = beta_p * R_v - alpha_p_v * p_v
    else:
        dp_v = 0

    dm_r = beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_r * m_r
    
    if (0 < p_tot < p_max):
        dp_r = beta_p * R_o * m_r/5e7  - alpha_p_r * p_r
        # 5 e 7 estimate of total amount of mRNA in cell
    else:
        dp_r = 0
  
    if (0 < R_v < R):
        dR_v = delta_R ** -1 * (dm_v * p_v/(1-p_v/p_max) + m_v * dp_v * p_max ** 2 / (p_max-p_v) ** 2  
    else:
        dR_v = 0
    
    dR_o = - dR_v
          
    return [dm_v, dp_v, dm_r, dp_r, dR_v, dR_o]

# list of initial conditions
inits = [m_v_init, p_v_init, m_r_init, p_r_init, R_v_init, R_o_init]

# time steps
t = np.linspace(0, 200, 1000)

# solving the ODE system
solution = ig.solve_ivp(model, [min(t),max(t)], inits, t_eval = t)
fig, axs = plt.subplot_mosaic([['mv', 'pv_vs_mv'],
                               ['pv', 'pv_vs_mv'],
                               ['mr', 'pr_vs_pv'],
                               ['pr', 'pr_vs_pv'],
                               ['fr', 'fr_vs_pv'],
                               ['gr', 'fr_vs_pv']])
axs['mv'].plot(solution.t, solution.y[0].T)
axs['mv'].set_title('mRNA of gene of value')
axs['pv'].plot(solution.t, solution.y[1].T)
axs['pv'].set_title('Protein of value')
axs['pv'].sharex(axs['mv'])
axs['mr'].plot(solution.t, solution.y[2].T)
axs['mr'].set_title('mRNA of regulator')
axs['mr'].sharex(axs['mv'])
axs['pr'].plot(solution.t, solution.y[3].T)
axs['pr'].set_title('Regulator protein')
axs['pr'].sharex(axs['mv'])
axs['fr'].plot(solution.t, solution.y[4].T)
axs['fr'].set_title('Value productive ribosomes')
axs['fr'].sharex(axs['mv'])
axs['gr'].plot(solution.t, solution.y[5].T)
axs['gr'].set_title('Other ribosomes')
axs['gr'].sharex(axs['mv'])
axs['gr'].set_xlabel('Time in hours')
axs['pv_vs_mv'].plot(solution.y[0].T, solution.y[1].T)
axs['pv_vs_mv'].set_title('state space prot val vs mRNA val')
axs['pr_vs_pv'].plot(solution.y[1].T, solution.y[3].T)
axs['pr_vs_pv'].set_title('state space reg prot vs prot val')
axs['fr_vs_pv'].plot(solution.y[1].T, solution.y[4].T)
axs['fr_vs_pv'].set_title('state space free ribosomes vs prot val')


plt.subplots_adjust(hspace=0.7)
plt.suptitle('Plots for delta_mv: ' + str(delta_m_v) + 'and delta_mr: ' + str(delta_m_r), fontsize=14)
plt.show()
