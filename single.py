import pyfasta


fasta = pyfasta.Fasta("/lustre/scratch118/malaria/team222/hh5/projects/analysis/assembly/anoph/Anoph_coluzzii_FALCON_Unzip_reseq_180829_primaryContigs.fasta")
keys_by_size = {key:len(fasta[key]) for key in fasta.keys()}
keys_by_size = sorted(keys_by_size.items(), key=lambda kv: -kv[1])
import subprocess
for (key,size) in keys_by_size:
    print key
    print size
    subprocess.call(["mummerplot","ref_qry.delta","-R","/lustre/scratch118/malaria/team222/hh5/ref/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa","-Q","/lustre/scratch118/malaria/team222/hh5/projects/analysis/assembly/anoph/Anoph_coluzzii_FALCON_Unzip_reseq_180829_primaryContigs.fasta","--layout" ,"--postscript","--filter","-q",key,"-p", key.replace("|","_")+str(size)])
