#script for making an anvio plot

#map reads from other metagenomic samples in the same diel cycle against the assembly of a metagenome using Bowtie2 and BamM
#make bowtie2 index
cat ($all_spades_assemblies) > all_cellular_contigs.fna
bowtie2-build all_cellular_contigs.fna cellular_contigs

for i in *fwd.fq.gz; do
bowtie2 \
--threads 16 \
-x cellular_contigs \
-1 $i \
-2 $(basename $i fwd.fq.gz).rev.fq.gz \
--no-unal \
-S $(basename $i fwd.fq.gz).sam 2>&1 | tee $(basename $i fwd.fq.gz).log;
done

#sort and index the bam files to order them
for i in *.sam; do
samtools view -F 4 -@ 16 -buSh $i | samtools sort -@ 16 - -o $(basename $i .sam).srt.bam;
samtools index -@ 16 $(basename $i .sam).srt.bam;
done

#filter mappings below >95% identity
for i in *.srt.bam; do
bamm filter -b $i --percentage_id 0.95 2>&1 | tee $(basename $i .bam).srt.fltr.log;
samtools sort -@ 16 $(basename $i .bam)_filtered.bam -o $(basename $i .bam).srt.fltr.bam;
samtools index $(basename $i .bam).srt.fltr.bam;
rm $(basename $i .bam)_filtered.bam ${i}*;
done

#generage the anvi'o database from the metagenomic assembly
anvi-gen-contigs-database -f all_cellular_contigs.fna -o contigs.db -n "AE1712 contigs database"
#run HMMs search for ribosomal RNA 
anvi-run-hmms -c contigs.db -T 16
#run COGs search to annotate protein families
anvi-run-ncbi-cogs -c contigs.db -T 16 --sensitive

#generate sequences for gene taxonomy
anvi-get-sequences-for-gene-calls -c contigs.db -o gene-calls.fa

#perform gene taxonomy on cellular contigs
KAIJUDB=$(kaiju_db_dir)

kaiju -t $KAIJUDB/nodes.dmp \
      -f $KAIJUDB/kaiju_db_nr.fmi \
      -i gene_calls.fa \
      -o gene_calls_nr.out \
      -z 16 \
      -v 2>&1 | tee -a $PROFILE_LOG

addTaxonNames -t $KAIJUDB/nodes.dmp \
              -n $KAIJUDB/names.dmp \
              -i gene_calls_nr.out \
              -o gene_calls_nr.names \
              -r superkingdom,phylum,class,order,family,genus,species

anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
-c contigs.db \
-p kaiju \
--just-do-it 2>&1 | tee -a $PROFILE_LOG
