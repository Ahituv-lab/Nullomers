import re,os,sys,glob
import mycustom2
from Bio.Seq import Seq 
from collections import Counter


def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

def fast_counter(SequencesL, kmer_lengthL): 
    Counts = Counter()
    for kmer_length in kmer_lengthL:
    	for seq in SequencesL:
        	for i in range(0, len(seq)+1-kmer_length):
        	    motif = seq[i:i+kmer_length]
        	    Counts[motif] += 1
    return Counts.most_common()


filed = "UP000005640_9606.fasta"
sequenceO = mycustom2.FastaFile(filed)
sequencesL = [ i.sequence.upper() for i in sequenceO ]

# Peptide lengths to scan
kmer_lengthL=range(1,8)
Counts=fast_counter(sequencesL,kmer_lengthL)

datafile=open("one_to_seven_amino_acids_UP000005640_9606","w")
for i in Counts:
	datafile.write(str(i[0])+'\t'+str(i[1])+'\n')
datafile.close()	
