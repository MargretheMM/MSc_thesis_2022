#!/usr/bin/env python3

import glob
import json

def score(prot, time):
    return float(prot) / float(time) ** 2

files = glob.glob("data_*.txt")

results = []
for filename in files:
    with open(filename) as handle:
        file_result ={}
        lines = handle.read().splitlines()
        _, file_result["p_v"], file_result["T2"] = lines[0].strip().split()
        file_result["score"] = score(file_result["p_v"], file_result["T2"])
        linenum = 2
        while linenum < len(lines):
            param_pairs = lines[linenum].split(":")
            file_result[param_pairs[0]] = float(param_pairs[1])
            linenum += 1
    if file_result not in results:
        results.append(file_result)

results.sort(key=lambda d: d["score"], reverse=True)

with open("top20scores.json", 'w') as handle:
    json.dump(results[:20],handle)
