#!/usr/bin/env python3

import sys
import glob


parameter_values = {}
for arg in sys.argv[1:]:
    k, v = arg.split("=")
    parameter_values[k] = v 
files = glob.glob("data_*.txt")

results = []
for filename in files:
    with open(filename) as handle:
        count = 0
        for line in handle.read().splitlines():
            if ":" not in line:
                continue
            name, value = line.split(": ")
            if name not in parameter_values:
                continue
            if value == parameter_values[name]:
                count += 1
    if count == len(parameter_values):
        results.append(filename)
    
for result in results:
    print(result)
