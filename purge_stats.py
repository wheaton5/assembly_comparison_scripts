

import pyfasta

fasta = pyfasta.Fasta("Anoph_coluzzii_FALCON_Unzip_reseq_180829_primaryContigsNOBAR.fasta")

num_haplotig = 0
num_repeat = 0
size_haplotig = 0
size_repeat = 0
size_curated = 0
num_curated = 0
num_junk = 0
size_junk = 0

with open("curated.reassignments.tsv") as curation:
    curation.readline()
    for line in curation:
        if line.startswith("#"):
            continue
        tokens = line.strip().split()
        curation_status = tokens[5]
        contig = tokens[0]
        if curation_status == "KEEP":
            num_curated += 1
            size_curated += len(fasta[contig])
        elif curation_status == "REPEAT":
            num_repeat += 1
            size_repeat += len(fasta[contig])
        elif curation_status == "HAPLOTIG":
            num_haplotig += 1
            size_haplotig += len(fasta[contig])
        elif curation_status == "JUNK":
            num_junk += 1
            size_junk += len(fasta[contig])
        else:
            print line
            assert(False) 

print "CURATED"
print str(num_curated) +" contigs"
print str(size_curated) +" bases"
print "REPEAT"
print str(num_repeat)+" contigs"
print str(size_repeat)+" bases"
print "HAPLOTIG"
print str(num_haplotig)+" contigs"
print str(size_haplotig)+" bases"
print "JUNK"
print str(num_junk)+" contigs"
print str(size_junk) + " bases"
