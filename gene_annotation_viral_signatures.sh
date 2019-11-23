prodigal -p meta -i <viral_input> -d <viral_output>

diamond blastx -d <diamond_nr_database> \
-q <prodigal_genecalls> \
-o <diamond_output> \
--outfmt 6 qseqid evalue bitscore stitle -k 1 \
--more-sensitive
