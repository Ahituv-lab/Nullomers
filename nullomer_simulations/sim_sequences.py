import re,os,sys,glob
import mycustom
import ushuffle
from Bio import SeqIO
import pdb

names=[]
sequencesL=[]

# path with the fasta file to be simulated 
filed = "/nfs/compgen-04/team218/ilias/nullomers_hg38_v2/hg38.fa"
sequenceO = mycustom.FastaFile(filed)
sequencesL = [ i.sequence.upper() for i in sequenceO ]
names = [ i.name.upper() for i in sequenceO ]

# Number of simulations
for k in range(1,101):
        datafile=open("sims_genome_dinucleotide/hg38_bootstrap_number_"+str(k)+"_controlling_dinucleotide_content.fasta","w")
        sequencesL_c=[]
        for index,i in enumerate(sequencesL):
                seq_random=ushuffle.shuffle(i,len(i), 2)
                datafile.write(">"+names[index]+'_control_bootstrap_'+str(k)+'\n')
                datafile.write(seq_random+'\n')
        datafile.close()
