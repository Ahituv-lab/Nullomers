Nullomers are extracted in three steps:

First we extract the kmer frequencies for kmers lengths provided with the command:
python extract_kmers.py

Secondly we produce all possible kmers of the same lengths with the command:
python possible_kmers.py

Third, we generate the final file with the frequency of all kmers and nullomers for those kmer lengths using the command:
python generate_all_counts.py


