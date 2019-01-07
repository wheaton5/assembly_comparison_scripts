

import pysam

#bam = pysam.AlignmentFile("all3_sorted.bam")
bam = pysam.AlignmentFile("all3_sorted.bam")
import numpy as np

#contig = "000066F_arrow_arrow"
#contig = "2L"
#start = 13411592
#end = 13427685
contig = "000057F_arrow_arrow"#"000022F_arrow_arrow"
#start = 13091592
#end = 14027685
#end = 29662443
length = bam.lengths[bam.gettid(contig)]
start = 0 
end = length
cov = np.zeros(end-start)
for read in bam.fetch(contig, start, end):
    #print "incrementing "+str(read.pos-start)+" to "+str(read.pos+read.alen-start)
    if read.mapq >= 60:
        cov[read.pos-start:read.pos+read.alen-start] += 1

with open("cov57F.csv",'w') as out:

    for i in range(start,end,5000):
        out.write(str(i)+","+str(np.mean(cov[i-start:min(end,i+5000-start)]))+"\n")
