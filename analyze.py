import pysam

bam = pysam.AlignmentFile("curated_alignments.bam")

supplemental_lengths = []
supplemental_distances = []
primary_lengths = []
unmapped_reads = 0
unmapped_bases = 0
with open("covered.bed",'w') as bed:
    for read in bam:
        chrom = read.reference_name
        pos = read.pos
        start = read.reference_start
        end = read.reference_end
        bed.write("\t".join([chrom,str(start),str(end)])+"\n")
        if read.is_supplementary or read.is_secondary:
            supplemental_lengths.append(end-start)
        else:
            primary_lengths.append(end-start)
        if read.is_unmapped:
            unmapped_reads += 1
            unmapped_bases += read.query_length

import numpy as np
print "primary"
print np.sum(primary_lengths)
print "supplemental"
print np.sum(supplemental_lengths)
print "unmapped"
print unmapped_reads
print unmapped_bases
