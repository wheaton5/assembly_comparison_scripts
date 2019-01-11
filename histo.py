import numpy as np

coverages = np.zeros(2000)

with open("coverage2.tsv") as cov:
    for line in cov:
        tok = line.strip().split("\t")
        x = int(tok[2])
        if x < 2000:
            coverages[x] += 1

with open("hist.csv",'w') as hist:
    for i,x in enumerate(coverages):
        hist.write(str(i)+","+str(x)+"\n")
