coverages = {}
import pyfasta
import numpy as np
fasta = pyfasta.Fasta("/lustre/scratch118/malaria/team222/hh5/ref/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa")
for chrom in fasta.keys():
    print chrom
    coverages[chrom] = np.zeros(len(fasta[chrom]))

with open("covered.bed") as bed:
    for line in bed:
        tokens = line.strip().split()
        chrom = tokens[0]
        start = int(tokens[1])
        end = int(tokens[2])
        coverages[chrom][start:end] += 1




with open("coverage_curated_regions.bed",'w') as bed:
    for chrom in coverages.keys():
        covs = coverages[chrom]
        start = 0
        stop = 0
        print chrom
        #for index, cov in enumerate(covs):
        while True:
            last_cov = covs[start]
            try:
                while last_cov == covs[stop]:
                    stop += 1
                bed.write("\t".join([chrom,str(start), str(stop), str(last_cov)])+"\n")
                last_cov = covs[stop]
                start = stop
            except:
                break

#counts = np.zeros(10)
#for chrom, coverage_list in coverages.iteritems():
#    for x in coverage_list:
#        if x < 10:
#            counts[x] += 1
totals = {'2L':48500000,'2R':60100000,'3L':40700000,'3R':52200000,'X':23400000,'UNKN':27200000}

total_counts = {}
for chrom, coverage_list in coverages.iteritems():
    if not chrom in totals:
        continue
    unique, counts = np.unique(coverage_list, return_counts=True)
    #print unique
    #print counts
    print chrom
    total = float(totals[chrom])#float(np.sum(counts))
    uncovered = 0
    covered = 0
    overcovered = 0
    for x,y in zip(unique, counts):
        if x == 0:
            print "uncovered "+str(y/total)
        elif x > 1:
            covered += y
            overcovered += y
        elif x == 1:
            covered += y
        total_counts.setdefault(x,0)
        total_counts[x] += y
    print "covered "+str(covered/total)
    print "overcovered "+str(overcovered/total)
print total_counts
covered = 0
overcovered = 0
for x, y in total_counts.iteritems():
    if x > 0:
        covered += y
    if x > 1:
        overcovered += y
print covered
print overcovered
