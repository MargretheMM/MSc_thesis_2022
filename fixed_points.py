#!/usr/bin/env python3

import numpy as np
import sympy as sp

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
# maximum translation in aas per hour per ribosome : aa/ribosome/sec * sec/hour = aa/hour
beta_p = 10 * 3600


# Set degradation rates and maximum production rates for mRNA and protein - mRNA: Curran et al 2013, Perez-Ortin et al 2007, protein Christiano 2020
# rates based on half-life in h
alpha_m_v = np.log(2)/0.33
alpha_p_v = np.log(2)/5


''' Presumed engineerable paramters '''
# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 1e6
K_p_r = 50

# mRNA production rates - mix of copy number, promotor strength etc - will be scaled with transcription rate
delta_m_v, delta_m_r = 1, 0.2

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
m_v, p_v, m_r, p_r, R_v, R_o = sp.symbols('m_v, p_v, m_r, p_r, R_v, R_o')

C = p_max*p_v/(p_max-p_v)

# equations

eq1 = sp.Eq(beta_m * delta_m_v *  control_negative(K_p_r, p_r)  - alpha_m_v * m_v,0)

eq2 = sp.Eq(beta_p * R_v - alpha_p_v * p_v,0)

eq3 = sp.Eq(beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_r * m_r,0)

eq4 = sp.Eq(beta_p * R_o * m_r/(5e7-m_v)  - alpha_p_r * p_r,0)

eq5 = sp.Eq(beta_p ** -1 * (eq1.lhs * C + m_v * eq2.lhs * p_max ** 2 / (p_max-p_v) ** 2 ),0)

eq6 = sp.Eq(R_o,R-R_v)

#solving
sol = sp.solve((eq1, eq2, eq3, eq4, eq5, eq6),(m_v, p_v, m_r, p_r, R_v, R_o))
print(sol)


