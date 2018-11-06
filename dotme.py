import pysam
import numpy as np

read_offsets = {}
read_positions = {}
read_chroms = {}
read_query_lengths = {}
contig_lengths = {}
assembly = "curated_assembly.fasta"
assembly_alignments = "curated_alignments.bam"
with open(assembly) as ass:
    current = 0
    name = None
    for line in ass:
        if line.startswith(">"):
            if current > 0:
                contig_lengths[name] = current
                current = 0
            name = line.strip().split(">")[1]
        else:
            current += len(line) - 1
    if current > 0:
        contig_lengths[name] = current

tid_to_chrom = {}
with open("dotplot.csv",'w') as dotplot:
    dotplot.write("chrom,ref_pos_start,ref_pos_end,read_pos_start,read_pos_end,read_name,direction,primary_chrom,primary_chrom_offset\n")
    read_directions = {}
    for bamname in [assembly_alignments]:
        bam = pysam.AlignmentFile(bamname)
        for read in bam:
            if read.is_unmapped:
                continue
            tid_to_chrom[read.tid] = bam.references[read.tid]
            read_directions.setdefault(read.qname,[0,0])
            read_offsets.setdefault(read.qname,{})
            read_chroms.setdefault(read.qname,{})
            query_length = read.query_alignment_end - read.query_alignment_start
            if read.cigartuples[0][0] == 5 or read.cigartuples[0][0] == 4:
                query_length += read.cigartuples[0][1]
            if read.cigartuples[-1][0] == 5 or read.cigartuples[-1][0] == 4:
                query_length += read.cigartuples[-1][1]
            read_query_lengths[read.qname] = query_length
            alen = read.alen
            read_offsets[read.qname].setdefault(read.tid,[])
            read_offsets[read.qname][read.tid].append((read.pos+alen/2.0, alen))
            read_chroms[read.qname].setdefault(read.tid, 0)
            read_chroms[read.qname][read.tid] += alen
            if read.is_reverse:
                read_directions[read.qname][1] += alen
            else:
                read_directions[read.qname][0] += alen
    major_chrom = {}
    with open("major_chrom.tsv",'w') as majorityreport:
        for (readname, tid_map) in read_chroms.iteritems():
            maxlen = -1
            secondbest = -1
            secondchrom = None
            max_chrom = None
            for (tid, alen) in tid_map.iteritems():
                if alen > maxlen:
                    secondbest = maxlen
                    secondchrom = max_chrom
                    maxlen = alen
                    max_chrom = tid
            if tid_to_chrom[max_chrom] == "UNKN" and not secondchrom == None:
                max_chrom = secondchrom
            major_chrom[readname] = max_chrom
            majorityreport.write("\t".join([readname,tid_to_chrom[max_chrom]])+"\n")
    average_position = {}
    denoms = {}
    for (readname, tid) in major_chrom.iteritems():
        denoms.setdefault(readname,0.0)
        #if tid_to_chrom[tid] == "UNKN":
        #    continue
        for (pos, length) in read_offsets[readname][tid]:
            denoms[readname] += float(length)
    for (readname, tid) in major_chrom.iteritems():
        average_position.setdefault(tid,{})
        average_position[tid].setdefault(readname,0)
        denoms.setdefault(readname,0.0)
        for (pos, length) in read_offsets[readname][tid]:
            average_position[tid][readname] += pos/denoms[readname]*length

    contig_order = {}
    for chrom in average_position.keys():
        contig_order[chrom] = [readname for (readname, avg_pos) in sorted(average_position[chrom].items(), key = lambda x: x[1])]
    read_offsets = {}
    for chrom in contig_order.keys():
        so_far = 0
        for readname in contig_order[chrom]:
            read_offsets[readname] = so_far 
            so_far += read_query_lengths[readname]


    for bamname in [assembly_alignments]:
        bam = pysam.AlignmentFile(bamname)

        for read in bam:
            if read.is_unmapped:
                continue
            overall_direction = True
            if read_directions[read.qname][1] > read_directions[read.qname][0]:
                overall_direction = False
            index = 0
            #if read.qname in read_offsets:
            #    if read_offsets[read.qname][1] > read.query_alignment_length:
            #        read_offsets[read.qname] = (read.pos-read.query_alignment_start, read.query_alignment_length)
            #else:
            #    read_offsets[read.qname] = (read.pos-read.query_alignment_start, read.query_alignment_length)
            offset = 0
            query_length = read.query_alignment_end - read.query_alignment_start
            if read.cigartuples[0][0] == 5 or read.cigartuples[0][0] == 4:
                query_length += read.cigartuples[0][1]
            if read.cigartuples[-1][0] == 5 or read.cigartuples[-1][0] == 4:
                query_length += read.cigartuples[-1][1]

            read_start = read.query_alignment_start 
            read_end = read.query_alignment_end
            if read.cigartuples[0][0] == 5:
                offset = read.cigartuples[0][1]
                read_start += offset
                read_end += offset
            if (read.is_reverse and overall_direction) or ((not read.is_reverse) and (not overall_direction)):
                tmp = read_start
                read_start = query_length - read_start
                read_end = query_length - read_end
            #if read.cigartuples[-1][0] == 5:
            #    read_end -= read.cigartuples[-1][1]
            ref_start = read.reference_start
            ref_end = read.reference_end
            read_dir = "forward"
            if (read.is_reverse and overall_direction) or ((not read.is_reverse) and (not overall_direction)):
                read_dir = "reverse"
            dotplot.write(",".join([bam.references[read.tid], str(ref_start), str(ref_end), str(read_start), str(read_end), read.qname, read_dir, bam.references[major_chrom[read.qname]], str(read_offsets[read.qname])])+"\n")


chrom_contig_lengths = {}
chroms = set()
for read, tid in major_chrom.iteritems():
    chrom = tid_to_chrom[tid]
    chroms.add(chrom)
    #print chrom
    chrom_contig_lengths.setdefault(chrom,[])
    chrom_contig_lengths[chrom].append(contig_lengths[read])
for chrom in chroms:
    lengths = sorted(chrom_contig_lengths[chrom], reverse=True)
    total_length = np.sum(lengths)
    print "total length of contigs for chrom "+chrom +" is  "+str(total_length)
    so_far = 0
    for contig_length in lengths:
        so_far += contig_length
        if so_far > 0.5*total_length:
            print "N50 for contigs for chrom "+chrom
            print contig_length
            break
