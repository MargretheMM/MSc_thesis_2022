#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt
import glob
import json
from numpy import mean

files = glob.glob("*.json")

shape = sys.argv[1]
var = ["m_v", "p_v", "m_r", "p_r", "T2"]
duration = 300

# get time series for plotting
for filename in files:
    if filename.split(":")[0] == shape:
        with open(filename) as handle:
            data = json.load(handle)

if not data:
    print("No relevant files found, please try again")
    sys.exit()

colours = ("#990000", "#2F3EEA",  "#1FD082",  "#FC7634",  "#030F4F")

data["time_hours"] = [x/3600 for x in data["time"]]

fig, axs = plt.subplot_mosaic([['m_v', 'm_v_zoom'],
                               ['p_v', 'p_v_zoom'],
                               ['m_r', 'm_r_zoom'],
                               ['p_r', 'p_r_zoom'],
                               ['T2', 'T2_zoom'],
                              ], figsize=(11,9))
index = 0
for v in var:
    axs[v].plot(data["time_hours"][:duration], data[v][:duration], color = colours[index])
    axs[v].set(xlabel='time (h)', ylabel = f"{v} in 10³ codons", title = f"{v} over time")
    axs[v].yaxis.set_major_formatter('{x:.3g}')
    if v == "T2":
        axs[v].set(ylabel = f"{v} in hours")
    axs[v].grid(color="lightgrey")
    stable_max = max(data[v][60:])
    stable_min = min(data[v][60:])
    stable_mean = mean(data[v][60:])
    fluctuations = stable_max - stable_min
    axs[v+"_zoom"].plot(data["time_hours"][:duration], data[v][:duration], color = colours[index])
    axs[v+"_zoom"].set(xlabel='time (h)', ylabel = f"{v} in 10³ codons", ylim=[stable_mean-fluctuations/1.5,stable_mean+fluctuations/1.5], title = f"{v} over time zoomed in")
    axs[v+"_zoom"].yaxis.set_major_formatter('{x:.3g}')
    if v == "T2":
        axs[v+"_zoom"].set(ylabel = f"{v} in hours")
    axs[v+"_zoom"].grid(color="lightgrey")
    index +=1

plt.subplots_adjust(hspace=0.7, wspace = 0.3)
fig.savefig(f"Plot_{shape}.png")
plt.show()
