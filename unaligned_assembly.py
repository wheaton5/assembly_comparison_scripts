import pysam


bam = pysam.AlignmentFile("curated_alignments.bam")


query_lengths = {}
query_covered = {} # map of query name to list of regions to be sorted and merged before taking stock

majority_chrom = {}
with open("major_chrom.tsv") as major:
    for line in major:
        tokens = line.strip().split()
        majority_chrom[tokens[0]] = tokens[1]

unaligned_sequence_by_chrom = {}

for read in bam:
    if read.qname == "000020F|arrow|arrow":
        print "known bacterial"
        continue
    query_length = read.query_alignment_end - read.query_alignment_start
    if read.cigartuples[0][0] == 5 or read.cigartuples[0][0] == 4:
        query_length += read.cigartuples[0][1]
    if read.cigartuples[-1][0] == 5 or read.cigartuples[-1][0] == 4:
        query_length += read.cigartuples[-1][1]
    if read.qname in query_lengths:
        assert(query_length == query_lengths[read.qname])
    query_lengths[read.qname] = query_length
    read_start = read.query_alignment_start
    read_end = read.query_alignment_end
    if read.cigartuples[0][0] == 5:
        offset = read.cigartuples[0][1]
        read_start += offset
        read_end += offset
    if read.is_reverse:
        tmp = read_end
        read_end = query_length - read_start
        read_start = query_length - tmp
    query_covered.setdefault(read.qname,[])
    query_covered[read.qname].append((read_start, read_end))

print "poo"
not_covered = 0
for query_name, regions in query_covered.iteritems():
    length = query_lengths[query_name]
    regs = sorted(regions)
    merged = []
    index = 0
    while index < len(regs):
        reg_start = regs[index][0]
        reg_end = regs[index][1]
        while index+1 < len(regs) and regs[index+1][0] < reg_end:
            reg_end = max(regs[index][1], regs[index+1][1])
            index += 1
        merged.append((reg_start, reg_end))
        index += 1
    covered = 0
    for region in merged:
        covered += region[1]-region[0]
    print merged
    print covered
    print length
    assert(length >= covered)
    print "query not covered for "+str(length-covered)+" bases"
    not_covered += length-covered
    chrom = majority_chrom[query_name]
    unaligned_sequence_by_chrom.setdefault(chrom,0)
    
    unaligned_sequence_by_chrom[chrom] += length-covered

print "by chrom"
print unaligned_sequence_by_chrom
     

print "total not covered "+str(not_covered)
