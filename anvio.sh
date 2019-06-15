#script for making an anvio plot

#map reads from other metagenomic samples in the same diel cycle against the assembly of a metagenome
bbwrap.sh ref=AE1712_3_80m_1k_reformatted.fasta build=3 append \
in=AE1712_1_interleaved_80m.qtrimmed.fq.gz,AE1712_1_interleaved_80m.merged.fq.gz,\
AE1712_5_interleaved_80m.qtrimmed.fq.gz,AE1712_5_interleaved_80m.merged.fq.gz,\
AE1712_3_interleaved_80m.qtrimmed.fq.gz,AE1712_3_interleaved_80m.merged.fq.gz,\
AE1712_7_interleaved_80m.qtrimmed.fq.gz,AE1712_7_interleaved_80m.merged.fq.gz,\
AE1712_9_interleaved_80m.qtrimmed.fq.gz,AE1712_9_interleaved_80m.merged.fq.gz,\
AE1712_11_interleaved_80m.qtrimmed.fq.gz,AE1712_11_interleaved_80m.merged.fq.gz \
out=AE1712_3_80m_1k_reformatted.sam

#sort and index the bam files to order them
samtools view -F 4 -buS AE1712_3_80m_1k_reformatted.sam | samtools sort - -o AE1712_3_80m_1k_reformatted.srt.bam
samtools index AE1712_3_80m_1k_reformatted.srt.bam

#generage the anvi'o database from the metagenomic assembly
anvi-gen-contigs-database -f AE1712_3_80m_1k_reformatted.fasta -o contigs_3.db -n "AE1712_3_80m contigs database"
#run HMMs search for ribosomal RNA 
anvi-run-hmms -c contigs_3.db -T 16
#run COGs search to annotate protein families
anvi-run-ncbi-cogs -c contigs_3.db -T 16 --sensitive

#generate sequences for gene taxonomy
anvi-get-sequences-for-gene-calls -c contigs_3.db -o gene-calls.fa

#use the Contig Annotation Tool to annotate the gene calls with taxonomy
CAT contigs \
-c gene_calls.fa \
-d ~/../bt273/BIOS-SCOPE/metag/ashley/tools/CAT/2018-10-16_CAT_database \
-t ~/../bt273/BIOS-SCOPE/metag/ashley/tools/CAT/2018-10-16_taxonomy \
-n 16 \
--out_prefix gene_calls_CAT

CAT add_names \
-i gene_calls_CAT.contig2classification.txt \
-o gene_calls_CAT.contig2classification.official_names.txt \
-t ~/../bt273/BIOS-SCOPE/metag/ashley/tools/CAT/2018-10-16_taxonomy \
--only_official

CAT summarise \
-c gene_calls.fa \
-i gene_calls_CAT.official_names.txt \
-o gene_calls_CAT.summary.txt

