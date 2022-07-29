#!/usr/bin/env python3

import glob


def score(prot, time):
    return float(prot) / float(time) ** 2


files = glob.glob("data_*.txt")

results = []
for filename in files:
    with open(filename) as handle:
        file_result =[]
        lines = handle.read().splitlines()
        _, proteins, time = lines[0].strip().split()
        file_result = [filename, proteins, time, str(score(proteins, time))]
        parameters = lines[2:]
        file_result.extend(parameters)
    results.append(file_result)

results.sort(key=lambda x: x[3], reverse=True)
for result in results[:30]:
    print(result)

with open("results.csv",'w') as handle:
    for result in results[:30]:
        handle.write(";".join(result)+"\n")
