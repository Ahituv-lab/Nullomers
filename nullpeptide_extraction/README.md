Nullpeptide are extracted in three steps:

First we extract the kmer frequencies for peptide kmers lengths provided with the command:
python extract_peptides.py

Secondly we produce all possible peptide kmers of the same lengths with the command:
python possible_kmers.py

Third, we generate the final file with the frequency of all peptide kmers and nullomers for those peptide kmer lengths using the command:
python generate_all_counts.py


