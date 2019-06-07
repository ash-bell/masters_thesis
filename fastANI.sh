#!/bin/bash
fastANI --rl ref_list.txt --ql ref_list.txt -o fastANI.output -t 16 --matrix
sed 's|/gpfs/ts0/projects/Research_Project-172179/metag/ashley/sags/reformatted/||g; s|/gpfs/ts0/projects/Research_Project-172179/metag/ashley/sags/other_sar11_genomes/use/||g; s/.reformatted.fasta//g; s/.fasta//g' fastANI.out > trimmed_fastANI.tsv
