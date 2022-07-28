#!/usr/bin/env python3

import sys
import scipy.integrate as ig
from matplotlib import pyplot as plt
import numpy as np




#    Model parameters - constants
#    Paramter values per cell, rates per second, concentration scale of 10^3 codons (per cell)
#        R_productive: Maximum number of productive ribosomes
#        p_max: Total amount of protein
#        m_max: Total amount of mRNA
#        beta_m: Maximum mRNA production rate per gene of strength 1
#        beta_p: Maximum protein production rate per ribosome
#        beta_R: Maximum ribosome production rate
#        alpha_m: General mRNA degradation rate
#        alpha_p: General protein degradation rate
#        alpha_R: Ribosome degradation rate
#        k: Maximum growth rate
#        K_m_v: Concentration constant for protein production based on mRNA concentration
#        m_v_init: Initial concentration of mRNA of value
#        p_v_init: Initial concentration of protein of value
#        m_r_init: Initial concentration of regulator mRNA
#        p_r_init: Initial concentration of regulator protein
#        R_init: Initial number of productive ribosomes

R_productive = 1.5e5
p_max = 3e7
m_max = 3e4
beta_m = 1/4
beta_p = 0.01
beta_R = 2000/60
alpha_m = np.log(2)/1200
alpha_p = np.log(2)/18000
alpha_R = np.log(2)/21600
k = np.log(2)/6000
K_m_v = 1e3
m_v_init = 0.05
p_v_init = 0.1
m_r_init = 0
p_r_init = 0
R_init = 1.2e5       
    
# Values to trial for engineerable parameters
K_p_r_trial = [0.1,1,10,50]
K_p_v_trial = [1e5,1e6,1e7]
delta_v_trial = [20,30]
delta_r_trial = [1e-5,1e-4,1e-3,1e-2]
multi_m_r_trial = [1,4,6]
multi_p_r_trial = [1,3,7,10]




def run_simulation(K_p_r, K_p_v, delta_r, delta_v, multi_m_r, multi_p_r): 
    '''Functions and model '''
    # Control functions -  Michaelis-menten-like
    def control_positive(K, substrate):
        assert substrate >= -1e-3, substrate
        return substrate/(K+substrate)

    def control_negative(K, substrate):
        assert substrate >= -1e-3, substrate
        return K/(K+substrate)
        
    alpha_m_r = alpha_m * multi_m_r
    alpha_p_r = alpha_p * multi_p_r


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
#        global step_number
#        step_number += 1
        m_v, p_v, m_r, p_r, R = init_deps

        # Ribosome fractions
        R_v = R * m_v/m_max
        R_o = R - R_v

        # protein level in model
        p_t = p_v + p_r + R * 12

        # current growth rate
        gamma = k * R_o / R_productive

        #check no negative values
        vals = np.append(init_deps,[R_v,R_o])
    #   steps.append((step_number, indep, vals))
    #    if step_number % 1000 == 0:
    #        pretty_print_labels()
    #        for step, timestamp, values in steps[-1:]:
    #            pretty_print_values(step, timestamp, values)
    #        print("previous rates")
    #        for values in rates[-1:]:
    #            pretty_print_values(None, None, values)        
        for thing in vals:
#            if thing < -1e-10:
#                pretty_print_labels()
#                for step, timestamp, values in steps[-5:]:
#                    pretty_print_values(step, timestamp, values)
#                print("previous rates")
#                for values in rates[-5:]:
#                    pretty_print_values(None, None, values)
            assert thing >= -1e-3, thing

# equations
        if (0 < m_v < m_max):
            dm_v = beta_m * delta_v * control_negative(K_p_r, p_r) - alpha_m * m_v
        elif m_v >= m_max:
            dm_v =  - alpha_m * m_v
        else:
            dm_v = beta_m * delta_v * control_negative(K_p_r, p_r)

        if (0 < p_t < p_max):
            dp_v = beta_p * m_v * 10 * control_positive(K_m_v,m_v) - (alpha_p + gamma) * p_v
        elif p_t >= p_max:
            dp_v =  -(alpha_p + gamma) * p_v
        else:
            dp_v = beta_p * m_v * 10 * control_positive(K_m_v,m_v)

        dm_r = beta_m * delta_r * control_positive(K_p_v,p_v) - alpha_m_r * m_r

        if (0 < p_t < p_max):
            dp_r = beta_p * R_o * m_r/(m_max - m_v)  - (alpha_p_r + gamma) * p_r
        elif p_t >= p_max:
            dp_r =  -(alpha_p_r + gamma) * p_r
        else:
            dp_r = beta_p * R_o * m_r/(m_max - m_v)

        if (0 < R < R_productive) and p_t < p_max:
            dR = beta_R * R_o / R_productive - (alpha_R + gamma) * R
        if R >= R_productive or p_t+dp_v >= p_max:
            dR =  -(alpha_R + gamma) * R
        else:
            dR = beta_R * R_o / R_productive

        res = [dm_v, dp_v, dm_r, dp_r, dR]
        return res

    # list of initial conditions
    inits = [m_v_init, p_v_init, m_r_init, p_r_init, R_init]

    # time steps
    t = np.linspace(0, 1e5, 10000)
    # 1.7e6 seconds is 20 days
    # 4.3e5 seconds is 5 days

    # solving the ODE system
    solution = ig.solve_ivp(model, [min(t),max(t)], inits, method='RK45', max_step = 0.3, t_eval = t)
    #Getting R_v and R_o from solution
    R_v_array = solution.y[4]*solution.y[0]/m_max
    R_o_array = solution.y[4] - R_v_array
    T2_array = np.log(2)/(k * R_o_array/R_productive) / 3600

    #plotting
    fig, axs = plt.subplot_mosaic([['mv', 'T2'],
                                   ['pv', 'textblock'],
                                   ['mr', 'textblock'],
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
    #axs['R'].set_ylim([1e5,1.7e5])
    axs['Rv'].plot(solution.t, R_v_array.T)
    axs['Rv'].set_title('Value productive ribosomes')
    axs['Rv'].sharex(axs['mv'])
    #axs['Rv'].set_ylim([0,2e5])
    axs['Ro'].plot(solution.t, R_o_array.T)
    axs['Ro'].set_title('Other ribosomes')
    axs['Ro'].sharex(axs['mv'])
    #axs['Ro'].set_ylim([0,2e5])
    axs['Ro'].set_xlabel(f'Time in seconds')
    axs['T2'].plot(solution.t, T2_array.T)
    axs['T2'].set_title('Doubling time in hours')
    axs['T2'].set_xlabel(f'Time in seconds')
    #axs['pv_vs_mv'].plot(solution.y[0].T, solution.y[1].T)
    #axs['pv_vs_mv'].set_title('state space prot val vs mRNA val')
    #axs['pr_vs_pv'].plot(solution.y[1].T, solution.y[3].T)
    #axs['pr_vs_pv'].set_title('state space reg prot vs prot val')
    #axs['R_vs_pv'].plot(solution.y[1].T, solution.y[4].T)
    #axs['R_vs_pv'].set_title('ribosomes vs prot val')

    constants = [
        [
            f'K_p_v: {K_p_v:.2g}',
            f'K_p_r: {K_p_r:.2g}',
            f'delta_v: {delta_v}',
            f'delta_r: {delta_r:.2g}',
            f'multi_m_r: {multi_m_r}',
            f'multi_p_r: {multi_p_r}',
        ],
        [
            f'beta_m: {beta_m:.3f}',
            f'beta_p: {beta_p:.3f}',
            f'beta_R: {beta_R:.3f}',
        ],
        [
            f'K_m_v: {K_m_v:.2g}',
            f'K_R: {K_R:.2g}',
            f'max_prots: {p_max:.2g}',
            f'max_prod_ribos: {R_productive:.2g}',
        ],
        [
            f'alpha_m: {alpha_m:.2g}',
            f'alpha_p: {alpha_p:.2g}',
            f'alpha_R: {alpha_R:.2g}',
        ]
    ]
    final_vals = [
                    f'Final mRNA of value: {solution.y[0].T[-1]:.3g}',
                    f'Final protein of value: {solution.y[1].T[-1]:.3g}',
                    f'Final regulator mRNA: {solution.y[2].T[-1]:.3g}',
                    f'Final regulator protein: {solution.y[3].T[-1]:.3g}',
                    f'Final number of ribosomes: {solution.y[4].T[-1]:.3g}',
                    f'Final Other ribosomes: {R_o_array.T[-1]:.3g}'
    ]
    constant_lines = "\n\n".join("\n".join(sublist) for sublist in constants)
    val_lines = "\n".join(final_vals)
    axs["textblock"].get_yaxis().set_visible(False)
    axs["textblock"].get_xaxis().set_visible(False)
    plt.figtext(.60,.15,constant_lines+"\n\n"+val_lines)
    plt.subplots_adjust(hspace=0.7)

    if __name__ == "__main__":
        if "--save" in sys.argv:
            import datetime
            name = datetime.datetime.now().isoformat().rsplit(".", 1)[0].replace(":", ".")
            plt.savefig(fname=f'plot_{name}.png')
            with open(f"data_{name}.txt",'w') as handle:
                handle.write("\t".join(["score", str(np.mean(solution.y[1].T[-100:])), str(np.mean(T2_array.T[-100:]))])+"\n")
                handle.write("values"+"\n")
                handle.write("\n".join(thing for thing in constants[0]))
                
        if "--noshow" not in sys.argv:
            plt.show()
    plt.close()             
        
# run model for combinations of parameters - only using degradation multipliers and K_p_r if also using regulation. 
for K_p_r in K_p_r_trial:
    for K_p_v in K_p_v_trial:
        for delta_v in delta_v_trial:
            for delta_r in delta_r_trial:
                for multi_m_r in multi_m_r_trial:
                    for multi_p_r in multi_p_r_trial:
                        try:
                            run_simulation(K_p_r,K_p_v, delta_r, delta_v, multi_m_r, multi_p_r)  
                        except AssertionError as err:
                            print(
                                err, '\nValues are:',
                                f'K_p_v: {K_p_v:g}', f'K_p_r: {K_p_r}',
                                f'delta_v: {delta_v}', f'delta_r: {delta_r}',
                                f'multi_m_r: {multi_m_r}', f'multi_p_r: {multi_p_r}'
                            )
                            continue

for K_p_v in K_p_v_trial:
    for delta_v in delta_v_trial:
        delta_r = 0
        multi_m_r = 0
        multi_p_r = 0
        K_p_r = 2
        try:
            run_simulation(K_p_r,K_p_v, delta_r, delta_v, multi_m_r, multi_p_r)  
        except AssertionError as err:
            print(
                err, '\nValues are:',
                f'K_p_v: {K_p_v:g}', f'K_p_r: {K_p_r}',
                f'delta_v: {delta_v}', f'delta_r: {delta_r}',
                f'multi_m_r: {multi_m_r}', f'multi_p_r: {multi_p_r}'
            )
            continue
                             
