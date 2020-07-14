# Nullomer finding code

This repository contains the code used in the manuscript [Absent from DNA and protein: genomic characterization of nullomers and nullpeptides across functional categories and evolution](https://www.biorxiv.org/content/10.1101/2020.03.02.972422v1) by Georgakopolous-Soares et al. The folders contain scripts for identifying and analyzing nullomers from multiple genomes. Please note that this is not intended to be a package for nullomer analysis but rather a place where the code that was used for the manuscript can be accessed and reviewed.

## nullomer_extractions



## nullomer_mutations

This directory contains [Julia](https://docs.julialang.org/en/v1/) files that are used to identify all possible nullomer creating single basepair deletion, insertion and substitution in the human genome.

* **findAlmostNullomers.jl** - this file scans the genome and for a given length it identifies all possible mutations creating kmers (does not have to be nullomers) that are created through a single basepair change. For the human genome it is recommended to use a high performance cluster for k>11. The flags -s and -e are meant to facilitate this by only considering a specific subset of the nullomers on the list.
* **findAlmostNullpeptides.jl** - similar to findAlmostNullomers.jl but operates on the proteome instead. Only reports substitutions. Again it is recommended to split up the process creating files that process 1,000 nullpeptides each using the -s and -e flags.
* **findNullomerMutations.jl** - this script takes as input a list of potential mutations (generated using findAlmostNullomers.jl) and a list of variants from gnomad.
* **gnomad_raw_extract.sh** - this script parses the (large) file from gnomad and extracts allele frequencies for all variants that involve a single bp substitution. This file is used by findNullomerMutations.jl. It will be useful to split the large file generated into smaller files, one per chromosome (using 23 for X and 24 for Y).
* **gnomad_afr_extract.sh** - this file extracts the variants that are population specific for afr. Minor modifications are required to extract variants specific for other populations.

Note that all of the above files will require editing to include suitable paths to genomes and other necessary resources.

### Identification of nucleotide substitutions resulting in nullpeptide creation

The findAlmostNullpeptides.jl script identifies all _amino acid_ substitutions that can result in the creation of a nullpeptide. However, due to constraints imposed by the genetic code, only a subset of these substitutions will be accessible through single basepair _nucleotide_ substitions. To identify this subset of mutations, additional processing is required:

* Run [annovar](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/) to identify the genomic regions that code for proteins
`perl /path/to/annovar/annotate_variation.pl -hgvs -buildver hg38 -dbtype wgEncodeGencodeBasicV33 /path/to/vcf/file/with/all/mutations/overlapping/coding/regions/all_cds_mutations.vcf /path/to/annovar/db/humandb/ -out /path/to/output/all_cds_mutations`
* Extract the non-synonymous changes `cut -f 2-8 /path/to/output/all_cds_mutations.exonic_variant_function | grep -e '^nonsyn' | cut -f 2-7 | awk {'printf ("%s\t%s\t%s\t%s\t%s\n", $2, $3, $1, $5, $6)'} > /path/to/output/all_cds_mutations_nonsyn.tsv`
* Process the output from findAlmostNullpeptides.jl to put into a smaller file. The shorter file is assumed to contain three columns with the transcript name (in Gencode/Ensembl format), the mutation (in AAposAA format, eg Q432W) and the index of the nullpeptide being created from the is mutation (refers to the index of the list used for findAlmostNullpeptides.jl) `grep ENST /path/to/file/nullpeptides5subspos.tsv | sed 's/\(ENST[[:digit:]]*\.[[:digit:]]*\)\t\([[:digit:]]*\)\t[A-Z]\([A-Z]\)[A-Z]:\([A-Z]\)\t[0-9]\t\([0-9]*\)/\1\t\3\2\4\t\5/' | sort -u > /path/to/file/nullpeptides5subspos_short.tsv`
* To identify the mutations that correspond to these mutations run a small shell script. As this takes a long time it is recommended to split into smaller jobs on a high performance cluster and belowed it is assumed that the coordinates for the nullpeptide generating mutations are stored in files where 1,000 nullpeptides at a time were processed

    n=`wc -l /path/to/nullpeptidesfile.txt`
    nfiles=$((n/1000 + 1))
    for ((i=1; i<=$nfiles; i++)) ; 
    do 
        <command for submitting job to scheduler> run_find_all_nullpeptide_mutations.sh $i ;
    done
* The results can then be merged as `cat /path/to/files/all_cds_nonsyn_nullpeptide5subs_* | sort -u > /path/to/files/all_cds_nonsyn_nullpeptide5subs_uniq.tsv`
* The total number of mutations is found as `wc -l /path/to/files/Gencode/all_cds_nonsyn_nullpeptide5subs_uniq.tsv`
* We can add back the id of the nullpeptide as well `grep -f /path/to/file/all_cds_nonsyn_one_transcript_nullpeptide5subs_uniq.tsv /path/to/file/nullpeptides5subspos_short.tsv > /path/to/file/all_cds_nonsyn_one_transcript_nullpeptide5subs_uniq_nullpepinds.tsv

## nullomer_permutations

## nullomer_simulation

## nullpeptide_extraction

