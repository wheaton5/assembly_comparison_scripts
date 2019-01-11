import pysam


#bam = pysam.AlignmentFile("m54119_180807_160930.subreads.bam",'rb',check_header=False, check_sq=False)
import gzip

import numpy as np
num_reads = 1000000
lengths = []#np.zeros(num_reads)
total = 0
with gzip.open("all.fq.gz") as fq:
    for index, line in enumerate(fq):
        if index/4 > num_reads:
            break
        if index % 4 == 1:
            lengths.append(len(line)-1)
            total += len(line)-1
lengths = sorted(lengths)
so_far = 0
done = True
for x in lengths:
    so_far += x
    if so_far >= 0.95*total:
        print x
        break



