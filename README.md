# Nullomer finding code

This repository contains the code used in the manuscript [Absent from DNA and protein: genomic characterization of nullomers and nullpeptides across functional categories and evolution](https://www.biorxiv.org/content/10.1101/2020.03.02.972422v1) by Georgakopolous-Soares et al. The folders contain scripts for identifying and analyzing nullomers from multiple genomes.

## nullomer_extractions



## nullomer_mutations

This directory contains [Julia](https://docs.julialang.org/en/v1/) files that are used to identify all possible nullomer creating single basepair deletion, insertion and substitution in the human genome.

* **findAlmostNullomers.jl** - this file scans the genome and for a given length it identifies all possible mutations creating kmers (does not have to be nullomers) that are created through a single basepair change. For the human genome it is recommended to use a high performance cluster for k>11. The flags -s and -e are meant to facilitate this by only considering a specific subset of the nullomers on the list.
* **findAlmostNullpeptides.jl** - similar to findAlmostNullomers.jl but operates on the proteome instead. Only reports substitutions.
* **findNullomerMutations.jl** - this script takes as input a list of potential mutations (generated using findAlmostNullomers.jl) and a list of variants from gnomad.
* **gnomad_raw_extract.sh** - this script parses the (large) file from gnomad and extracts allele frequencies for all variants that involve a single bp substitution. This file is used by findNullomerMutations.jl
* **gnomad_afr_extract.sh** - this file extracts the variants that are population specific for afr. Minor modifications are required to extract variants specific for other populations.

Note that all of the above files will require editing to include suitable paths to genomes and other necessary resources.

## nullomer_permutations

## nullomer_simulation

## nullpeptide_extraction

