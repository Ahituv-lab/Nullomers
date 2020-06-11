import re,os,sys,glob
import mycustom 
from Bio.Seq import Seq 
from collections import Counter
import pickle


def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def fast_counter(SequencesL,SequencesL_rev_compl, kmer_lengthL): 
    Counts = Counter()
    for kmer_length in kmer_lengthL:
    	for seq in SequencesL:
        	for i in range(0, len(seq)+1-kmer_length):
        	    motif = seq[i:i+kmer_length]
        	    Counts[motif] += 1
        for seq2 in SequencesL_rev_compl:
                for im in range(0, len(seq2)+1-kmer_length):
                    motif = seq2[im:im+kmer_length]
                    Counts[motif] += 1
    return Counts.most_common()


filed = "hg38.fasta"
sequenceO = mycustom.FastaFile(filed)
sequencesL = [ i.sequence.upper() for i in sequenceO ]
sequencesL_rev_compl = []
for i in sequencesL:
        seq = Seq(i)
        sequencesL_rev_compl+=[str(seq.reverse_complement())]

# kmer lengths to scan - here e scan 1bp to 15bp long kmers
kmer_lengthL=range(1,16)
Counts=fast_counter(sequencesL,sequencesL_rev_compl,kmer_lengthL)
	
datafile=open("hg38_data","w")
for i in Counts:
	datafile.write(str(i[0])+'\t'+str(i[1])+'\n')
datafile.close()	
