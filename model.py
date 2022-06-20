#!/usr/bin/env python3

import sys
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np


''' Model parameters - biologically based '''

# maximum number of ribosomes available for protein production - currently assume 95% of total - Metzl-Raz 2017 and von der Haar 2008
R_max= 2e5 * 0.95


# total number of RNA polymerases in cell - Borggrefe 2001
total_pol = 30000

# Maximum number of protein bound amino acids in cell - futcher et al 1999
p_max = 3e10
# Total codons of mRNA in cell - estimate from Siewiak 2010 - take high esitmate
m_tot = 3e7

# Maximum transcription and translation rates
# max transcription in codons per hour per cell : codons/sec * sec/hour *number of polymerases = codons/hour - Yu and Nielsen 2019
beta_m = 50/3 * 3600 * total_pol
# maximum translation in aas per hour per ribosome : aa/ribosome/sec * sec/hour = aa/hour
beta_p = 10 * 3600


# max rate of production of new ribosomes : max rate per minute * min/hour =  max s per hour - Warner 1999
beta_R = 2000 * 60

# Set degradation rates and maximum production rates for mRNA and protein - mRNA: Curran et al 2013, Perez-Ortin et al 2007, protein Christiano 2020
# rates based on half-life in h
alpha_m_v = np.log(2)/0.33
alpha_p_v = np.log(2)/5
alpha_R = np.log(2)/6

#Maximum growth rate (based on doubling time of 80 minutes)
k=np.log(2)*3/4

''' Presumed engineerable paramters '''
# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 1e6
K_p_r = 50
K_m_v = 1e2
K_R = 0.1 * R_max

# mRNA production rates - mix of copy number, promotor strength etc - will be scaled with transcription rate
delta_m_v, delta_m_r = 0.5, 0

# regulator degradation rates
alpha_m_r = np.log(2)/0.33
alpha_p_r = np.log(2)/1

''' Initial conditions'''
# initial amounts of each product (mRNA and protein from each gene)
m_v_init = 5
p_v_init = 100
m_r_init = 0
p_r_init = 0
R_init = R_max * 0.99
R_v_init = 5
R_o_init = R_init - R_v_init


'''Functions and model '''
# Control functions -  Michaelis-menten-like
def control_positive(K, substrate):
    return substrate/(K+substrate)

def control_negative(K, substrate):
    return K/(K+substrate)

# tracking if growth stops (R_o hits zero)
def no_growth(t,y): return y[6]-1
# stop simulation if this happens
#no_growth.terminal = True

# write up equations
def model(indep: float, init_deps):
    '''
    Defining the set of differential equations to be solved
    Input: independent variable, initial values of dependent variables (list)
    Output: list of equations in the system 
    '''
    m_v, p_v, m_r, p_r, R, R_v, R_o = init_deps
    
    # protein level in model
    p_tot = p_v + p_r + R * 12e3
    
    # current growth rate
    gamma = k * R_o / R
    
    # equations
    if (0 < m_v < m_tot):
        dm_v = beta_m * delta_m_v *  control_negative(K_p_r, p_r)  - alpha_m_v * m_v
    else:
        dm_v = 0
            
    if (0 < p_tot < p_max):
        dp_v = beta_p * control_positive(K_R,R) * control_positive(K_m_v,m_v) - (alpha_p_v + gamma) * p_v
    else:
        dp_v = 0

    dm_r = beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_r * m_r
    
    if (0 < p_tot < p_max):
        dp_r = beta_p * R_o * m_r/(m_tot - m_v)  - (alpha_p_r + gamma) * p_r

    else:
        dp_r = 0
  
    if (0 < R < R_max):
        dR = beta_R * R_o / R_max - (alpha_R + gamma) * R
    else:
        dR = 0
    
    dR_v = m_tot ** -1 * (dR * m_v + R * dm_v)
    
    dR_o = dR - dR_v

    return [dm_v, dp_v, dm_r, dp_r, dR, dR_v, dR_o]

# list of initial conditions
inits = [m_v_init, p_v_init, m_r_init, p_r_init, R_init, R_v_init, R_o_init]

# time steps
t = np.linspace(0, 200, 100)

# solving the ODE system
solution = ig.solve_ivp(model, [min(t),max(t)], inits, t_eval = t, events = no_growth)
#, events=no_growth)
fig, axs = plt.subplot_mosaic([['mv', 'pv_vs_mv'],
                               ['pv', 'pv_vs_mv'],
                               ['mr', 'pr_vs_pv'],
                               ['pr', 'pr_vs_pv'],
                               ['R', 'R_vs_pv'],
                               ['Rv', 'R_vs_pv'],
                               ['Ro', 'R_vs_pv']])
axs['mv'].plot(solution.t, solution.y[0].T)
axs['mv'].set_title('mRNA of gene of value')
axs['pv'].plot(solution.t, solution.y[1].T)
axs['pv'].set_title('Protein of value')
axs['pv'].sharex(axs['mv'])
axs['pv'].set_ylim([0,3e6])
axs['mr'].plot(solution.t, solution.y[2].T)
axs['mr'].set_title('mRNA of regulator')
axs['mr'].sharex(axs['mv'])
#axs['mr'].set_ylim([0,3e7])
axs['pr'].plot(solution.t, solution.y[3].T)
axs['pr'].set_title('Regulator protein')
axs['pr'].sharex(axs['mv'])
#axs['pr'].set_ylim([0,3e10])
axs['R'].plot(solution.t, solution.y[4].T)
axs['R'].set_title('Ribosomes')
axs['R'].sharex(axs['mv'])
axs['R'].set_ylim([0,2e5])
axs['Rv'].plot(solution.t, solution.y[5].T)
axs['Rv'].set_title('Value productive ribosomes')
axs['Rv'].sharex(axs['mv'])
axs['Rv'].set_ylim([0,2e5])
axs['Ro'].plot(solution.t, solution.y[6].T)
axs['Ro'].set_title('Other ribosomes')
axs['Ro'].sharex(axs['mv'])
axs['Ro'].set_ylim([0,2e5])
axs['Ro'].set_xlabel('Time in hours')
axs['pv_vs_mv'].plot(solution.y[0].T, solution.y[1].T)
axs['pv_vs_mv'].set_title('state space prot val vs mRNA val')
axs['pr_vs_pv'].plot(solution.y[1].T, solution.y[3].T)
axs['pr_vs_pv'].set_title('state space reg prot vs prot val')
axs['R_vs_pv'].plot(solution.y[1].T, solution.y[4].T)
axs['R_vs_pv'].set_title('ribosomes vs prot val')


plt.subplots_adjust(hspace=0.7)
plt.suptitle('Plots for delta_mv: ' + str(delta_m_v) + 'and delta_mr: ' + str(delta_m_r), fontsize=14)
plt.show()
