#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt
import glob
import json

        
files = glob.glob("*.json")

plot_data = {}
times = None

# get time series for plotting
rank = 1
for filename in files:
    with open(filename) as handle:
        data = json.load(handle)
        if times is None:
            times = data["time"]
        if rank < 11:
            plot_data[str(rank)] = data["p_v"]
        rank += 1
            

if not plot_data:
    print("No relevant files found, please try again")
    sys.exit()
    

colours = ("#990000", "#2F3EEA", "#1FD082", "#FC7634","#030F4F", "#008835", "#79238E","#99000080", "#2F3EEA80", "#1FD08280", "#FC763480","#030F4F80", "#00883580", "#79238E80")

fig, axs = plt.subplots(2,figsize=(11,10))
time_in_hours = [x/3600 for x in times]

stable_min = 3e7
stable_max = 0

for key in plot_data.keys():
    axs[0].plot(time_in_hours,plot_data[key], color = colours[int(key)], label = f"rank: {key}")
    axs[1].plot(time_in_hours,plot_data[key],  color = colours[int(key)], label = f"rank: {key}")
    dummy = plot_data[key]
    long_term = dummy[-100:]
    if stable_min > min(long_term):
        stable_min = min(long_term) 
    if stable_max < max(long_term):
        stable_max = max(long_term)
        
fluctuations = stable_max-stable_min

axs[0].set(xlabel='time (h)', ylabel = f"p_v in 10³ codons", title = f"p_v over time of top 10 scoring parameter sets")
axs[0].yaxis.set_major_formatter('{x:.3g}')       
axs[0].grid(color="lightgrey")
axs[0].legend(loc='lower right')
axs[1].set(xlabel='time (h)', ylabel = f"p_v in 10³ codons", ylim=[stable_min-fluctuations/10,stable_max+fluctuations/10], title = f"p_v over time of top 10 scoring parameter sets zoomed in")
axs[1].yaxis.set_major_formatter('{x:.5g}')
axs[1].grid(color="lightgrey")
axs[1].legend(loc='lower left')
plt.subplots_adjust(wspace = 0.1)
fig.savefig(f"Plot_thirdit.png")
plt.show()                    
