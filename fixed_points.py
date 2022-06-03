#!/usr/bin/env python3

import numpy as np
import sympy as sp

''' Model parameters - biologically based '''

# total number of ribosomes - von der Haar 2008
total_s = 2e5
# maximum number of ribosomes available for protein production - currently assume 95% of total - Metzl-Raz 2017
s_max = total_s*0.95
# half full stock
s_half = s_max / 2

# total number of RNA polymerases in cell - Borggrefe 2001
total_pol = 30000


# Maximum transcription and translation rates
# max transcription in codons per hour per cell : codons/sec * sec/hour *number of polymerases = codons/hour - Yu and Nielsen 2019
beta_m = 50/3 * 3600 * total_pol
# maximum translation in aas per hour per cell : aa/ribosome/sec * sec/hour * number of ribsomes which can be diverted = aa/hour
delta_s = 10 * 3600
beta_p = delta_s * s_max

# max rate of production of new ribosomes : max rate per minute * min/hour =  max s per hour - Warner 1999
beta_s = 2000 * 60

# Set degradation rates and maximum production rates for mRNA and protein and ribosomes - mRNA: Curran et al 2013, Perez-Ortin et al 2007, protein Christiano 2020
# rates based on half-life in h
alpha_m_v = np.log(2)/0.33
alpha_p_v = np.log(2)/5
alpha_s = np.log(2)/6


# proportionality constants
# protein level influence on ribosome production
delta_p = 1e-8
# ribosome number influence on growth rate
delta_gamma = 0.1

''' Presumed engineerable paramters '''
# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 1e6
K_p_r = 50

# mRNA production rates - mix of copy number, promotor strength etc - will be scaled with transcription rate
delta_m_v, delta_m_r = 1, 1

# regulator degradation rates
alpha_m_r = np.log(2)/0.33
alpha_p_r = np.log(2)/1


'''Functions and model '''
# Control functions. Currently Michaelis-menten-like (Hill functions, so far n = 1)
def control_positive(K, substrate):
    return substrate/(K+substrate)

def control_negative(K, substrate):
    return K/(K+substrate)


# defining symbols
m_v, p_v, m_r, p_r, s, gamma = sp.symbols('m_v, p_v, m_r, p_r, s, gamma')

# equations

eq1 = sp.Eq(beta_m * delta_m_v *  control_negative(K_p_r, p_r)  - alpha_m_v * m_v,0)

eq2 = sp.Eq(beta_p * control_positive(s_half, s) * m_v - (alpha_p_v + gamma) * p_v,0)

eq3 = sp.Eq(beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_r * m_r,0)

eq4 = sp.Eq(beta_p * control_positive(s_half, s) * m_r - (alpha_p_r + gamma) * p_r,0)

eq5 = sp.Eq(beta_s * gamma - delta_s ** -1 * m_v - delta_p * (p_v + p_r) - (alpha_s + gamma) * s,0)

eq6 = sp.Eq(delta_gamma * (s/s_half - 1),0)

#solving
sol = sp.solve((eq1, eq2, eq3, eq4, eq5, eq6),(m_v, p_v, m_r, p_r, s, gamma))
print(sol)


