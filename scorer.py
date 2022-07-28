#!/usr/bin/env python3

import glob


def score(prot, time):
    return float(prot) / float(time)


files = glob.glob("data_*.txt")

values_by_filename ={}
results = []
for filename in files:
    with open(filename) as handle:
        lines = handle.read().splitlines()
        _, proteins, time = lines[0].strip().split()
        parameters = lines[2:]
    results.append([filename, proteins, time, score(proteins, time)])
    values_by_filename[filename] = parameters
    
results.sort(key=lambda x: x[3], reverse=True)
for result in results[:10]:
    print(result, values_by_filename[result[0]])
