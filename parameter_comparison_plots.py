#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt
import glob
import json

        
files = glob.glob("*.json")

par = sys.argv[1]
var = "p_v"

plot_data = {}
times = None

# get time series for plotting
for filename in files:
    if filename.split(":")[0] == par:
        with open(filename) as handle:
            data = json.load(handle)
            if times is None:
                times = data["time"]
            param_value = data[par]
            time_series = data[var]
            plot_data[param_value] = time_series

if not plot_data:
    print("No relevant files found, please try again")
    sys.exit()

colours = ("#990000", "#2F3EEA", "#1FD082", "#FC7634","#030F4F", "#008835", "#79238E","#99000080", "#2F3EEA80", "#1FD08280", "#FC763480","#030F4F80", "#00883580", "#79238E80")

fig, ax = plt.subplots(figsize=(10,7.5))
index = 0
time_in_hours = [x/3600 for x in times]
for key in sorted(plot_data.keys()):
        ax.plot(time_in_hours[:300],plot_data[key][:300],color = colours[index], label = f"{par}: {key}")
        index += 1

ax.set(xlabel='time (h)', ylabel = "protein of value amount in 10³ codons")
ax.grid(color="lightgrey")
ax.legend(loc='lower right')

fig.savefig(f"Plot_{var}_vs_{par}_{round(time_in_hours[299])}.png")
plt.show()                    
