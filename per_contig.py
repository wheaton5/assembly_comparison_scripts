contig = "000059F|arrow|arrow"


import pyfasta
fasta = pyfasta.Fasta("Anoph_coluzzii_FALCON_Unzip_reseq_180829_primaryContigs.fasta")

with open(contig.replace("|","_")+".fa",'w') as out:
    out.write(">"+contig+"\n")
    out.write(fasta[contig][:])
    out.write("\n")

