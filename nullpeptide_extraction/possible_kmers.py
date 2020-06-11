import re,os,sys,itertools

def kmer_generate(lengthsL):
	total_kmers=[]
	AminoAcidsL=["G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T"]
	for k in lengthsL:
		total_kmers+=[''.join(p) for p in itertools.product(AminoAcidsL, repeat=k)]
	return total_kmers

Peptide_KmersL=kmer_generate(range(1,8))
		
datafile=open("peptide_kmers_1_7nt","w")
for i in Peptide_KmersL:
        datafile.write(str(i)+'\n')
datafile.close()  
