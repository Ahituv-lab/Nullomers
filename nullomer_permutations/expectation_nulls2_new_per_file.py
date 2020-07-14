import re,os,sys,glob,math
import numpy as np

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def kmer_reader(s):
        return s.split()[0],s.split()[-1]

def reader(path):
	OccsNull={}
	Occs={}
	with open(path) as f1: 
                for line in f1: 
                        motif,occs= kmer_reader(line)
			if int(occs)==0:
                        	OccsNull[motif]=int(occs)
			Occs[motif]=int(occs)
	return Occs,OccsNull

def motif_find(mot):
	letters=["A","G","C","T"]
	motifs_to_score=[]
	for k in range(len(mot)):
		for letter in letters:
			motif = mot[:k]+letter+mot[k+1:]
			motifs_to_score.append(motif)
	return list(set(motifs_to_score)-set(motif))

files = ["All_kmer_words_1_15nt_occurrences_CCDS_GENCODE_V28_hg38_15bp","All_kmer_words_1_15nt_occurrences_GENCODE_V28_exonic_hg38_15bp","All_kmer_words_1_15nt_occurrences_genic_GENCODE_V28_hg38_15bp","All_kmer_words_1_15nt_occurrences_5UTR_exons_GENCODE_V28_hg38_total_15bp","All_kmer_words_1_15nt_occurrences_3UTR_exons_GENCODE_V28_hg38_total_15bp","All_kmer_words_1_15nt_occurrences_intronic_GENCODE_V28_hg38_15bp","three_to_fifteenmers_GENCODE_V28_genes_hg38_genic_promoters_minus_2500_plus_500_total_15bp","All_kmer_words_1_15nt_occurrences_genome_hg38_15bp"]
picked = int(sys.argv[1])
files=[files[picked]]
for filed in files:
	datafile=open(filed+".perm","w")
	names=["coding","exonic","genic","5UTR","3UTR","intronic","promoters","genome"]
	Dic,DicNull=reader(filed)

	for k in DicNull.keys():
		motifs_expected=motif_find(k)
		expected = np.median([Dic[j] for j in motifs_expected])/float(len(motifs_expected))
		datafile.write(k+'\t'+str(expected)+'\n')
	datafile.close()

