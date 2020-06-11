import re,os,sys,glob
import mycustom 
from Bio.Seq import Seq 
from collections import Counter
import pickle
import io


def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

DataL=reader('one_to_seven_amino_acids_UP000005640_9606')
DataD = {}
for i in DataL:
        if "N" not in i[0]:
                DataD[i[0]]=int(i[1])


DataL=reader("peptide_kmers_1_7nt")
Data2D={}
for k in DataL:
        Data2D[k[0]]=0

DataTD= { k: DataD.get(k, 0) + Data2D.get(k, 0) for k in set(DataD) | set(Data2D) }

datafile=open("one_to_seven_amino_acids_UP000005640_9606_total","w")
for i in DataTD.items():
        datafile.write(i[0]+'\t'+str(i[1])+'\n')
datafile.close()
