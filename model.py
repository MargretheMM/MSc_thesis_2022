#!/usr/bin/env python3

import sys
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np

# initial amounts of each product (mRNA and protein from each gene) and resource level (maybe number of ribosomes?)
m_v_init = 10
p_v_init = 100
m_r_init = 2
p_r_init = 1

# total number of ribosomes
total_ribos = 190000
# maximum number of ribosomes available for protein production - currently assume 80% of total
max_avail_ribos = total_ribos*0.8

#initally, most ribosomes are available
free_ribos_init = max_avail_ribos*0.9

# Maximum transcription and translation rates
# transcription in codons per hour per cell : codons/sec * sec/hour * number of polymerases that can work at same time on same DNA = codons/hour
transcription_rate = 40/3 * 3600 * 3
# maximum translation in aas per hour per cell : aa/ribosome/sec * sec/hour * number of ribsomes which can be diverted = aa/hour
translation_rate_per_ribosome = 10 * 3600
max_translation_rate = translation_rate_per_ribosome * max_avail_ribos

# max rate of production of new ribosomes : max rate per minute * min/hour =  max ribos per hour
max_ribos_manufacturing_rate = 2000 * 60

# Set degradation rates and maximum production rates for mRNA and protein and ribosome (still needs checking)
# rates based on half-life in h, ^-1
deg_rate_m = 2
deg_rate_p = 0.2
deg_rate_ribos = 0.15

# mRNA production rates - mix of copy number, promotor strength etc - will be scaled with transcription rate
production_m_v, production_m_r = 1, 1

# initial cell growth rate - in generations per hour - use as fitness parameter, if falls too low -> collapse
dilution_rate_init = 0.5
max_dilution_rate = 1

# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 1000000
K_p_r = 50



# Control functions. Currently Michaelis-menten-like (Hill functions, so far n = 1)
def control_positive(K, substrate):
    return substrate/(K+substrate)

def control_negative(K, substrate):
    return K/(K+substrate)



# write up equations
def model(indep, init_deps):
    '''
    Defining the set of differential equations to be solved
    Input: independent variable, initial values of dependent variables (list)
    Output: list of equations in the system 
    '''
    m_v, p_v, m_r, p_r, ribos, dilution = init_deps

    # equations
    dm_v = transcription_rate * production_m_v *  control_negative(K_p_r, p_r)  - deg_rate_m * m_v
    
    dp_v = max_translation_rate * control_positive(K_stock, ribos) * m_v - (deg_rate_p + dilution) * p_v
    
    dm_r = transcription_rate * production_m_r * control_positive(K_p_v,p_v) - deg_rate_m * m_r
    
    dp_r = max_translation_rate * control_positive(K_stock, ribos) * m_r - (deg_rate_p + dilution) * p_r
    
    dribos = max_ribos_manufacturing_rate * dilution - translation_rate_per_ribosome ** -1 * m_v - (deg_rate_ribos + dilution) * ribos
    
    if dilution > 0 and dilution < max_dilution_rate:
        ddilution = 10**-6 * (ribos/(max_avail_ribos/2) - 1)
    else:
        ddilution = 0
          
    return [dm_v, dp_v, dm_r, dp_r, dribos, ddilution]

# list of initial conditions
inits = [m_v_init, p_v_init, m_r_init, p_r_init, free_ribos_init, dilution_rate_init]

# time steps
t = np.linspace(0, 200, 1000)

# solving the ODE system
solution = ig.solve_ivp(model, [min(t),max(t)], inits, t_eval = t)
plt.figure(1)
plt.subplot(611)
plt.plot(solution.t, solution.y[1].T)
plt.title('POI over time')
plt.subplot(612)
plt.plot(solution.t, solution.y[0].T)
plt.title('mRNA of gene 1')
plt.subplot(613)
plt.plot(solution.t, solution.y[2].T)
plt.title('mRNA of gene 2')
plt.subplot(614)
plt.plot(solution.t, solution.y[3].T)
plt.title('regulator protein')
plt.subplot(615)
plt.plot(solution.t, solution.y[4].T)
plt.title('free ribosomes')
plt.subplot(616)
plt.plot(solution.t, solution.y[5].T)
plt.title('growth rate)

#plt.legend(['m1', 'm2', 'p2', 'fitness index'])
plt.show()
