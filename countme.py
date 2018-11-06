total = 0
chrom_lengths = {}
with open("../../../../ref/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa") as ass:
    #with open("Anoph_coluzzii_FALCON_Unzip_reseq_180829_primaryContigs.fasta") as ass:
    #with open("curated.artefacts.fasta") as ass:
    #with open("curated_assembly.fasta") as ass:
    last = None
    count = 0
    for line in ass:
        if line.startswith(">"):
            if count > 0:
                chrom_lengths[last] = count
                total += count
            last = line.strip().split(">")[1]
            count = 0
        else:
            count += len(line.strip().replace("N","").replace("n",""))
    if count > 0:
        chrom_lengths[last] = count
        total += count
print chrom_lengths
print total
            
