#!/usr/bin/env python3

import sys
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np

# initial amounts of each product (mRNA and protein from each gene)
m_v_init = 10
p_v_init = 1000
m_r_init = 1
p_r_init = 1

# total number of ribosomes
total_s = 190000
# maximum number of ribosomes available for protein production - currently assume 95% of total
s_max = total_s*0.95

#initally, most ribosomes are available
s_init = s_max*0.7

# Maximum transcription and translation rates
# transcription in codons per hour per cell : codons/sec * sec/hour = codons/hour
beta_m = 40/3 * 3600
# maximum translation in aas per hour per cell : aa/ribosome/sec * sec/hour * number of ribsomes which can be diverted = aa/hour
delta_s = 10 * 3600
beta_p = delta_s * s_max

# max rate of production of new ribosomes : max rate per minute * min/hour =  max s per hour
beta_s = 2000 * 60

# Set degradation rates and maximum production rates for mRNA and protein and ribosome (still needs checking)
# rates based on half-life in h, ^-1
alpha_m = 2
alpha_p = 0.2
alpha_s = 0.15

# mRNA production rates - mix of copy number, promotor strength etc - will be scaled with transcription rate
delta_m_v, delta_m_r = 1000000, 1

# initial cell growth rate - in generations per hour - use as fitness parameter, if falls too low -> collapse
gamma_init = 0.5
max_gamma = 1

# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 1000000
K_p_r = 50
s_half = s_max / 2

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
    m_v, p_v, m_r, p_r, s, gamma = init_deps

    # equations
    dm_v = beta_m * delta_m_v *  control_negative(K_p_r, p_r)  - alpha_m * m_v
    
    dp_v = beta_p * control_positive(s_half, s) * m_v - (alpha_p + gamma) * p_v
    
    dm_r = beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m * m_r
    
    dp_r = beta_p * control_positive(s_half, s) * m_r - (alpha_p + gamma) * p_r
    
    if (0 < s < s_max):
        ds = beta_s * gamma - delta_s ** -1 * m_v - (alpha_s + gamma) * s
    else:
        ds = 0
    
    if (0 < gamma < max_gamma):
        dgamma = 0.1* (s/s_half - 1)
    else:
        dgamma = 0
          
    return [dm_v, dp_v, dm_r, dp_r, ds, dgamma]

# list of initial conditions
inits = [m_v_init, p_v_init, m_r_init, p_r_init, s_init, gamma_init]

# time steps
t = np.linspace(0, 200, 1000)

# solving the ODE system
solution = ig.solve_ivp(model, [min(t),max(t)], inits, t_eval = t)
plt.figure(1)
plt.subplot(611)
plt.plot(solution.t, solution.y[0].T)
plt.title('mRNA of gene of value')
plt.subplot(612)
plt.plot(solution.t, solution.y[1].T)
plt.title('Protein of value')
plt.subplot(613)
plt.plot(solution.t, solution.y[2].T)
plt.title('mRNA of regulator')
plt.subplot(614)
plt.plot(solution.t, solution.y[3].T)
plt.title('Regulator protein')
plt.subplot(615)
plt.plot(solution.t, solution.y[4].T)
plt.title('Free ribosomes')
plt.subplot(616)
plt.plot(solution.t, solution.y[5].T)
plt.title('Growth rate')
plt.ylim(0, 1)

#plt.legend(['m1', 'm2', 'p2', 'fitness index'])
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.7)
plt.suptitle('Plots for delta_mv: ' + str(delta_m_v) + 'and delta_mr: ' + str(delta_m_r), fontsize=14)
plt.show()
