import re,os,sys,itertools

def kmer_generate(lengthsL):
	total_kmers=[]
	bases=['A','T','G','C']
	for k in lengthsL:
		total_kmers+=[''.join(p) for p in itertools.product(bases, repeat=k)]
	return total_kmers

KmersL_1_15=kmer_generate(range(1,16))
		
datafile=open("kmers_1_15nt","w")
for i in KmersL_1_15:
        datafile.write(str(i)+'\n')
datafile.close()  
