#!/bin/bash
wrapper_phage_contigs_sorter_iPlant.pl \
	-f all_sags.fasta \
	--db 2 \
	--ncpu 32 \
	--data-dir ../../tools/virsorter/virsorter-data/ \
	--diamond \
	--no_c \
	2>&1 | tee virsorter.log
