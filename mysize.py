import pysam

files = ["m54119_180806_194558.subreads.bam",
"m54119_180807_160930.subreads.bam",
"m54119_180808_103633.subreads.bam"]
total_size = 0
for fn in files:
    bam = pysam.AlignmentFile(fn, 'rb',  check_header=False, check_sq=False)
    for index,read in enumerate(bam):
        total_size += len(read.seq)
        if index % 1000000 == 0:
            print total_size
            print index
            print fn
