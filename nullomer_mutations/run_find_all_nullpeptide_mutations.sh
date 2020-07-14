#!/bin/bash
i=$1
s=$(((i-1)*1000+1))
e=$((i*1000))
n=`wc -l /path/to/nullpeptides/nullpeptides5subspos_short_${s}_${e}.tsv | cut -f 1 --delim=' '` #this file is assumed to contain three columns with the transcript name (in Gencode/Ensembl format), the mutation (in AAposAA format, eg Q432W) and the index of the nullpeptide being created from the is mutation (refers to the index of the list used for findAlmostNullpeptides.jl). It can be created by straightforward use of awk from the output of findAlmostNullpeptides.jl
sl=1
el=1000
mutfilename="/path/to/coding/mutations/file/all_cds_mutations_nonsyn.tsv"
savefilename="/path/to/directory/all_cds_nonsyn_nullpeptide5subs_"${s}"_"${e}".tsv"

sed -n "${sl},${el}p;${el}q" /path/to/nullpeptides/nullpeptides5subspos_short.tsv | cut -f 1,2 | sed 's/\t/.*/' | grep -f - ${mutfilename} > ${savefilename}
nm=$((n/1000 + 1))
for ((j=1; j<=nm; j++))
do
    sl=$(((j-1)*1000+1))
    el=$((j*1000))
    sed -n "${sl},${el}p;${el}q" /path/to/nullpeptides/nullpeptides5subspos_short_${s}_${e}.tsv | cut -f 1,2 | sed 's/\t/.*/' | grep -f - ${mutfilename} >> ${savefilename}
done
