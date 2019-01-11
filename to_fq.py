import pysam


bam = pysam.AlignmentFile("m54119_180808_103633.subreads.bam",'rb', check_header=False, check_sq=False)
with open("m54119_180808_103633.fq",'w') as fq:
    for read in bam:
        fq.write("@"+read.qname+"\n")
        fq.write(read.seq+"\n")
        fq.write("+\n")
        fq.write(read.qual+"\n")
