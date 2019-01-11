sequence = []
import numpy as np
with open("../Anoph_coluzzii_FALCON_Unzip_reseq_180829_curated_primary.clean.noBact.wMT.split15F.sort.fasta") as anoph:
    current = 0
    for line in anoph:
	if line.startswith(">"):
            if current > 0:
                sequence.append(current)
                current = 0
	    continue
	current += len(line.strip())
    if current > 0:
        sequence.append(current)

sequence = sorted(sequence)
print sequence
with open("blah.csv",'w') as blah:
    for x in sequence:
        blah.write(str(x)+"\n")
print np.sum(sequence)
total = float(np.sum(sequence))
sequence_so_far = 0.0
for x in sequence:
    sequence_so_far += x
    if sequence_so_far >= total*0.5:
        print x
        break
