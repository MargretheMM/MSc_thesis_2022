#!/usr/bin/env python3

import sys
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np



''' Model parameters - biologically based '''
# time base - the base timestep for the various constants in seconds
t_base = 1

# maximum number of ribosomes available for protein production - currently assume 75% of total - Metzl-Raz 2017 and von der Haar 2008
R_max = 2e5
R_productive = R_max * 0.75


# total number of RNA polymerases in cell - Borggrefe 2001 - can't acutally use all polymerases on same gene, serial, not parallel.
# total_pol = 30000

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
beta_R = 2000 / 60 * t_base

# Set degradation rates and maximum production rates for mRNA and protein - mRNA: Curran et al 2013, Perez-Ortin et al 2007, protein Christiano 2020
# rates based on half-life in minutes
alpha_m_v = np.log(2)/1200 * t_base
alpha_p_v = np.log(2)/18000 * t_base
alpha_R = np.log(2)/21600 * t_base

#Maximum growth rate (based on doubling time of 80 minutes)
k=np.log(2)/4800 * t_base

''' Presumed engineerable paramters '''
# contcentration constants - association/disassociation in one - need some ideas for this
K_p_v = 5e6
K_p_r = 2
K_m_v = 5e3
K_R = 0.3 * R_productive

# mRNA production rates - mix of copy number, promotor strength etc - scales transcription rate
delta_m_v, delta_m_r = 10, 0

# regulator degradation rate multipliers
multi_m_r = 1
multi_p_r = 5
alpha_m_r = alpha_m_v * multi_m_r
alpha_p_r = alpha_p_v * multi_p_r

''' Initial conditions'''
# initial amounts of each product (mRNA and protein from each gene)
m_v_init = 0.05
p_v_init = 0.1
m_r_init = 0
p_r_init = 0
R_init = R_productive * 0.80


'''Functions and model '''
# Control functions -  Michaelis-menten-like
def control_positive(K, substrate):
    assert substrate >= -1e-5, substrate
    return substrate/(K+substrate)

def control_negative(K, substrate):
    assert substrate >= -1e-5, substrate
    return K/(K+substrate)

# tracking if growth stops (R hits zero)
def no_growth(t,y): return y[4]-1
# stop simulation if this happens
no_growth.terminal = True

# other slow growth check
def too_slow(rate):
    assert rate > np.log(2)/36000 * t_base, "culture grows too slowly"
    
steps = []
rates = []

step_number = 0

def pretty_print_labels():
    parts = "step, time, R, R_v, R_o".split(", ")
    parts = ["%9s" % key for key in parts]
    print("".join(parts))

def pretty_print_values(step, timestamp, values):
    ts = ("%7.2f" % timestamp) if timestamp is not None else " " * 9
    vals = ["%1.3g" % value for value in values]
    print(("%9d" % step) if step else " "*7, ts, "".join("%9s" % value for value in vals))


# write up equations
def model(indep: float, init_deps):
    '''
    Defining the set of differential equations to be solved
    Input: independent variable, initial values of dependent variables (list)
    Output: list of equations in the system
    '''
    global step_number
    step_number += 1
    m_v, p_v, m_r, p_r, R = init_deps

    # Ribosome fractions
    R_v = R * m_v/m_tot
    R_o = R - R_v

    # protein level in model
    p_tot = p_v + p_r + R * 12

    # current growth rate
    gamma = k * R_o / R_productive
    too_slow(gamma)
    #check no negative values
    vals = np.append(init_deps[-1],[R_v,R_o])
    steps.append((step_number, indep, vals))
#    if step_number % 1000 == 0:
#        pretty_print_labels()
#        for step, timestamp, values in steps[-1:]:
#            pretty_print_values(step, timestamp, values)
#        print("previous rates")
#        for values in rates[-1:]:
#            pretty_print_values(None, None, values)        
    for thing in vals:
        if thing < -1e-10:
            pretty_print_labels()
            for step, timestamp, values in steps[-5:]:
                pretty_print_values(step, timestamp, values)
            print("previous rates")
            for values in rates[-5:]:
                pretty_print_values(None, None, values)
            assert thing >= -1e-10, thing
#    if R > R_max:
#        pretty_print_labels()
#        for step, timestamp, values in steps[-5:]:
#            pretty_print_values(step, timestamp, values)
#        print("previous rates")
#        for values in rates[-5:]:
#            pretty_print_values(None, None, values)
#            print(alpha_R, gamma, "->", -(alpha_R + gamma) * R)
#        raise

    # equations
    if (0 < m_v < m_tot):
        dm_v = beta_m * delta_m_v * control_negative(K_p_r, p_r) - alpha_m_v * m_v
    elif m_v >= m_tot:
        dm_v =  - alpha_m_v * m_v
    else:
        dm_v = 0

    if (0 < p_tot < p_max):
        dp_v = beta_p * R_v * control_positive(K_R,R) * control_positive(K_m_v,m_v) - (alpha_p_v + gamma) * p_v
    elif p_tot >= p_max:
        dp_v =  -(alpha_p_v + gamma) * p_v
    else:
        dp_v = 0

    dm_r = beta_m * delta_m_r * control_positive(K_p_v,p_v) - alpha_m_r * m_r

    if (0 < p_tot < p_max):
        dp_r = beta_p * R_o * m_r/(m_tot - m_v)  - (alpha_p_r + gamma) * p_r
    elif p_tot >= p_max:
        dp_r =  -(alpha_p_r + gamma) * p_r
    else:
        dp_r = 0

    if (0 < R < R_productive) and p_tot < p_max:
        dR = beta_R * R_o / R_productive - (alpha_R + gamma) * R
    if R >= R_productive or p_tot+dp_v >= p_max:
        dR =  -(alpha_R + gamma) * R
    else:
        dR = 0

    res = [dm_v, dp_v, dm_r, dp_r, dR]
    rates.append(res)
    return res

# list of initial conditions
inits = [m_v_init, p_v_init, m_r_init, p_r_init, R_init]

# time steps
t = np.linspace(0, 4.3e5/t_base, 1000)
# 1.7e7 seconds is 20 days
# 4.3e5 seconds is 5 days

# solving the ODE system
solution = ig.solve_ivp(model, [min(t),max(t)], inits, method='BDF', max_step = 0.5, t_eval = t, events = no_growth)
#Getting R_v and R_o from solution
R_v_array = solution.y[4]*solution.y[0]/m_tot
R_o_array = solution.y[4] - R_v_array
#plotting
fig, axs = plt.subplot_mosaic([['mv', 'phase'],
                               ['pv', 'phase'],
                               ['mr', 'phase'],
                               ['pr', 'textblock'],
                               ['R', 'textblock'],
                               ['Rv', 'textblock'],
                               ['Ro', 'textblock'],
                              ],figsize=(11,10))
axs['mv'].plot(solution.t, solution.y[0].T)
axs['mv'].set_title('mRNA of gene of value')
axs['pv'].plot(solution.t, solution.y[1].T)
axs['pv'].set_title('Protein of value')
axs['pv'].sharex(axs['mv'])
#axs['pv'].set_ylim([0,3e6])
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
#axs['R'].set_ylim([0,2e5])
axs['Rv'].plot(solution.t, R_v_array.T)
axs['Rv'].set_title('Value productive ribosomes')
axs['Rv'].sharex(axs['mv'])
#axs['Rv'].set_ylim([0,2e5])
axs['Ro'].plot(solution.t, R_o_array.T)
axs['Ro'].set_title('Other ribosomes')
axs['Ro'].sharex(axs['mv'])
#axs['Ro'].set_ylim([0,2e5])
axs['Ro'].set_xlabel(f'Number of time steps of length {t_base} seconds')
#axs['pv_vs_mv'].plot(solution.y[0].T, solution.y[1].T)
#axs['pv_vs_mv'].set_title('state space prot val vs mRNA val')
#axs['pr_vs_pv'].plot(solution.y[1].T, solution.y[3].T)
#axs['pr_vs_pv'].set_title('state space reg prot vs prot val')
#axs['R_vs_pv'].plot(solution.y[1].T, solution.y[4].T)
#axs['R_vs_pv'].set_title('ribosomes vs prot val')

constants = [
    [
        f'Time unit: {t_base}s',
        f'max_prots: {p_max:.2g}',
        f'max_prod_ribos: {R_productive:.2g}',
    ],
    [
        f'beta_m: {beta_m:.3f}',
        f'beta_p: {beta_p:.3f}',
        f'beta_R: {beta_R:.3f}',
        f'delta_mv: {delta_m_v}',
        f'delta_mr: {delta_m_r}',
    ],
    [
        f'K_m_v: {K_m_v}',
        f'K_R: {K_R:g}',
        f'K_p_v: {K_p_v:g}',
        f'K_p_r: {K_p_r}',
    ],
    [
        f'alpha_m_v: {alpha_m_v:.2g}',
        f'alpha_p_v: {alpha_p_v:.2g}',
        f'alpha_R: {alpha_R:.2g}',
        f'multi_m_r: {multi_m_r}',
        f'multi_p_r: {multi_p_r}',
    ]
]
constant_lines = "\n\n".join("\n".join(sublist) for sublist in constants)
axs["textblock"].get_yaxis().set_visible(False)
axs["textblock"].get_xaxis().set_visible(False)
plt.figtext(.66,.15,constant_lines)
plt.subplots_adjust(hspace=0.7)

if __name__ == "__main__":
    if "--save" in sys.argv:
        import datetime
        name = datetime.datetime.now().isoformat().rsplit(".", 1)[0].replace(":", ".")
        plt.savefig(fname=f'plot_{name}.png')
    if "--noshow" not in sys.argv:
        plt.show()
