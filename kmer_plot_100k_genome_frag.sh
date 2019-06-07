#!/bin/bash

for i in *.fasta; do
pyfasta split -n 15 -k 100000 $i;
done

cat *100kmer > 100k_genomes.fasta

python ~/tools/marine-phage-paper-scripts/kmer_freq.py 100k_genomes.fasta > 100k_genomes_kmer_count.fasta
python ~/tools/marine-phage-paper-scripts/run_umap.py -l 1 -d 0 -n 2 100k_genomes_kmer_count.fasta
python ~/tools/marine-phage-paper-scripts/run_hdbscan.py -p 100k_genomes -c 75 100k_genomes.umap.tsv
python ~/tools/marine-phage-paper-scripts/run_hdbscan.py -p 100k_genomes 100k_genomes.hdbscan.tsv

