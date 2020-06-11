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

DataL=reader('hg38_data')
DataD = {}
for i in DataL:
        if "N" not in i[0]:
                DataD[i[0]]=int(i[1])


DataL=reader("kmers_1_15nt")
Data2D={}
for k in DataL:
        Data2D[k[0]]=0

DataTD= { k: DataD.get(k, 0) + Data2D.get(k, 0) for k in set(DataD) | set(Data2D) }

datafile=open("hg38_data_total","w")
for i in DataTD.items():
        datafile.write(i[0]+'\t'+str(i[1])+'\n')
datafile.close()
