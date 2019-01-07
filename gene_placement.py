import pysam

bam = pysam.AlignmentFile("gene_assembly_alignments.bam")

gene_to_alignments = {}

for read in bam:
    name = read.qname
    gene_to_alignments.setdefault(name, [])
    gene_to_alignments[name].append(read)


def get_aligned_bases(read):
    query_length = read.query_alignment_end - read.query_alignment_start
    if read.cigartuples == None:
        print query_length
        print read.qname
        print read.tid
        print read.is_unmapped
        return 0
    if read.cigartuples[0][0] == 5 or read.cigartuples[0][0] == 4:
        query_length -= read.cigartuples[0][1]
    if read.cigartuples[-1][0] == 5 or read.cigartuples[-1][0] == 4:
        query_length -= read.cigartuples[-1][1]
    return query_length

contig_primary_chroms = {}
with open("dotplot.csv") as dot:
    dot.readline()
    for line in dot:
        tokens = line.strip().split(",")
        contig_primary_chroms[tokens[5]] = tokens[7]
        

with open("gene_placement.csv",'w') as out:
    out.write(",".join(["gene,chrom,contig\n"]))
    for gene, reads in gene_to_alignments.iteritems():
        extents = []
        aligned_bases = 0
        query_length = 0
        primary_contig = None
        contig_aligned_bases = {}
        for read in reads:
            #if read.is_unmapped:
            #   out.write(",".join([gene.split("_")[1],"None","None","False"])+"\n")
            if not read.is_supplementary:
                non_ns = 0
                seq = str(read.query_sequence)
                for c in seq:
                    if not (c == "N" or c == "n"):
                        non_ns +=1
            query_length = non_ns#max(query_length, len(read.seq))
            contig_aligned_bases.setdefault(read.tid,0)
            contig_aligned_bases[read.tid] += get_aligned_bases(read)
        primary_tid = sorted(contig_aligned_bases.items(), key=lambda kv: kv[1])[0][0]
        primary_contig = bam.references[primary_tid]
        print primary_contig
        if gene.split("_")[0] =="UNKN":
            if contig_aligned_bases[primary_tid]/float(query_length) > 0.95:
                out.write(",".join([gene.split("_")[1],contig_primary_chroms[primary_contig],primary_contig.split("_")[0]])+"\n")
            else:
                out.write(",".join([gene.split("_")[1],"not_placed","not_placed"])+"\n")
            
            #aligned_bases += query_end-query_start
            
