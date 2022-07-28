#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt
import glob
import json

        
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
                              ], figsize=(11,10))
index = 0
for v in var:
    axs[v].plot(data["time_hours"][:duration], data[v][:duration], color = colours[index])
    axs[v].set(xlabel='time (h)', title = f"{v} over time")
    axs[v].grid(color="lightgrey")
    index += 1
    fluctuations = max(data[v][50:]) - min(data[v][50:])
    axs[v+"_zoom"].plot(data["time_hours"][:duration], data[v][:duration], color = colours[index])
    axs[v+"_zoom"].set(xlabel='time (h)', ylim=[min(data[v][50:])-fluctuations,max(data[v][50:])+fluctuations], title = f"{v} over time zoomed in")
    axs[v+"_zoom"].grid(color="lightgrey")
    index +=1

plt.subplots_adjust(hspace=0.7)
fig.savefig(f"Plot_{shape}.png")
plt.show()                    
