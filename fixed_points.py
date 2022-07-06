#!/usr/bin/env python3

import numpy as np
import sympy as sp

''' Model parameters - biologically based '''
# time base - the base timestep for the various constants in seconds
t_base = 1

# maximum number of ribosomes available for protein production - currently assume 75% of total - Metzl-Raz 2017 and von der Haar 2008
R_max = 2e5
R_productive = R_max * 0.75

# Maximum number of protein bound amino acids in cell - futcher et al 1999
p_max = 3e7
# Total codons of mRNA in cell - estimate from Siewiak 2010 - take high esitmate
m_tot = 3e4

# Maximum transcription and translation rates
# max transcription in codons per polymerase * 10 (assumed number of polymerases which can work on same bit of DNA) = codons/sec Yu and Nielsen 2019
beta_m = 50/300 * t_base
# alternative way is calculating max rate similarly to max ribosome production rate, but using half-life: 3e7/2 codons should be produced in 15 minutes, which is about 17 000 codons per second - but this is cell wide -
# maximum translation in aas per ribosome : aa/ribosome/sec
beta_p = 0.01 * t_base

# max rate of production of new ribosomes : max rate per second  =  max ribos per minute / seconds per minute - Warner 1999
beta_R = 2100 / 60 * t_base

# Set degradation rates and maximum production rates for mRNA and protein - mRNA: Curran et al 2013, Perez-Ortin et al 2007, protein Christiano 2020
# rates based on half-life in minutes
alpha_m_v = np.log(2)/1200 * t_base
alpha_p_v = np.log(2)/18000 * t_base
alpha_R = np.log(2)/21600 * t_base

#Maximum growth rate (based on doubling time of 80 minutes)
k=np.log(2)/4800 * t_base

''' Presumed engineerable paramters '''
# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 1e6
K_p_r = 2
K_m_v = 1e2
K_R = 0.5 * R_productive

# mRNA production rates - mix of copy number, promotor strength etc - scales transcription rate
delta_m_v, delta_m_r = 5, 0.002

# regulator degradation rate multipliers
multi_m_r = 1
multi_p_r = 5

'''Functions and model '''
# Control functions. Currently Michaelis-menten-like (Hill functions, so far n = 1)
def control_positive(K, substrate):
    return substrate/(K+substrate)

def control_negative(K, substrate):
    return K/(K+substrate)


# defining symbols
m_v, p_v, m_r, p_r, R_v, R_o = sp.symbols('m_v, p_v, m_r, p_r, R_v, R_o')

# equations

eq1 = sp.Eq(beta_m * delta_m_v *  control_negative(K_p_r, p_r)  - alpha_m_v * m_v,0)

eq2 = sp.Eq(beta_p * R_v * m_v * control_positive(K_R,R) * control_positive(K_m_v,m_v) - (alpha_p_v + gamma) * p_v,0)

eq3 = sp.Eq(beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_r * m_r,0)

eq4 = sp.Eq(beta_p * m_r * R_o * m_r/(m_tot - m_v)  - (alpha_p_r + gamma) * p_r,0)

eq5 = sp.Eq(dR = beta_R * R_o / R_productive - (alpha_R + gamma) * R,0)

#solving
sol = sp.solve((eq1, eq2, eq3, eq4, eq5, eq6),(m_v, p_v, m_r, p_r, R_v, R_o))
print(sol)


