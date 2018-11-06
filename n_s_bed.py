import pyfasta

fasta = pyfasta.Fasta("/lustre/scratch118/malaria/team222/hh5/ref/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa")

with open("n_s_bed.bed",'w') as bed:
    for key in fasta.keys():
        chrom = fasta[key]
        print key
        start_index = 0
        index = 0
        #for index,base in enumerate(chrom):
        while True:
            try:                
                start_index = index
                while chrom[index] == 'n' or chrom[index] == 'N':
                    index += 1
                if not index == start_index:
                    bed.write(key+"\t"+str(start_index)+"\t"+str(index)+"\n")
                start_index = index
                while not (chrom[index] == 'n' or chrom[index] == 'N'):
                    index += 1
            except:
                break
