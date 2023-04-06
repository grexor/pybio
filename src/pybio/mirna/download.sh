rm *.fasta
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz -O mirna.fasta.gz
gunzip mirna.fasta.gz
python split_by_species.py
