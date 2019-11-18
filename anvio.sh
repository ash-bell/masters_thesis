#script for making an anvio

#rename all contigs with simple names
for i in *200m*.fasta; do
anvi-script-reformat-fasta \
$i \
-o $(basename $i .fasta).fna \
-l 2500 \
--simplify-names;
done

#generage the anvi'o database from the metagenomic assembly
anvi-gen-contigs-database -f all_cellular_contigs.fna -o contigs.db -n "Cellular contigs database"
#run HMMs search for ribosomal RNA 
anvi-run-hmms -c contigs.db -T 16
#run COGs search to annotate protein families
anvi-run-ncbi-cogs -c contigs.db -T 16 --sensitive

for i in *.bam; do
anvi-profile \
--input-file $i \
--contigs-db contigs.db \
--output-dir $(basename $i .bam)_PROFILE \
--sample-name $(basename $i .bam) \
--num-threads 16;
done

anvi-merge *_PROFILE/PROFILE.db -o SAMPLES-MERGED -c contigs.db
#generate sequences for gene taxonomy
anvi-get-sequences-for-gene-calls -c contigs.db -o gene-calls.fa

#make bowtie2 index
#cat all the 200m contigs together
cat *200m*.fna > all_cellular_contigs.fna
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

#perform gene taxonomy on cellular contigs
KAIJUDB=$(kaiju_db_dir)

kaiju -t $KAIJUDB/nodes.dmp \
      -f $KAIJUDB/kaiju_db_nr.fmi \
      -i gene_calls.fa \
      -o gene_calls_nr.out \
      -z 16 \
      -v 2>&1 | tee -a kaiju.db

addTaxonNames -t $KAIJUDB/nodes.dmp \
              -n $KAIJUDB/names.dmp \
              -i gene_calls_nr.out \
              -o gene_calls_nr.names \
              -r superkingdom,phylum,class,order,family,genus,species

anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
-c contigs.db \
-p kaiju \
--just-do-it 2>&1 | tee -a $PROFILE_LOG

#Get BinSanity bin identities per contig
get-ids -f . -l gene_calls.fa -o ids.txt -x 1
Binsanity-profile -i gene_calls.fa -s . --ids ids.txt -c gene_calls_coverage -T 16
Binsanity-wf -f . -l gene_calls.fa -c gene_calls_coverage.cov.x100.lognorm --threads 16 -o BinSanityWF_output

#Get MetaBAT2 bin identities per contig using recomended settings from their manual
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
metabat \
-i gene_calls.fa \
-a depth.txt \
-o metabat2 \
--maxP 90 \
--maxEdges 500 \
--minS 75 \
--noAdd \
--numThreads 16 | tee metabat2.log

#Get Maxbin bin identities per contig
for i in *.bam; do
pileup.sh in=$i out=$(basename $i .bam).pileup;
cut -f1,2 $(basename $i .bam).pileup > $(basename $i .bam)_cov_for_maxbin.tsv;
done

run_MaxBin.pl \
-thread 16 \
-contig gene_calls.fa \
-out maxbin \
-abund <1>_cov_for_maxbin.tsv \
-abund2 <2>_cov_for_maxbin.tsv \
-abund3 <3>_cov_for_maxbin.tsv \
-abund4 <4>_cov_for_maxbin.tsv \
-abund5 <5>_cov_for_maxbin.tsv \
-abund6 <6>_cov_for_maxbin.tsv | tee -a maxbin.log

#run in python to tabulate outputs of binning software into format anvi'o can read
```
#!/usr/bin/env python
from Bio import SeqIO #from biopython import the following packages SeqIO - helps read and write sequence data
import glob #helps with pathing on linux machines
import os #read and write files
import pandas as pd #additional data management options 



def process_file(results_list, output_name): #make a function called process_file with inputs to be results_list and output_name
    df = pd.DataFrame(results_list, columns=['contig', 'bin']) #object df is a dataframe with input results_list and contains columns defined as contigs and bins
    df.to_csv(output_name, sep='\t', index=False) #function to turn the df into a csv file seperated by tab without an indexing data

results = [] #make the results object 

#maxbin
file_list = glob.glob('./maxbin/*.fasta') #file_list defined as anything in maxbin filepath with fasta in it using glob package for pathing
for f in file_list: #for every item in the file list do
    bin_name = os.path.basename(f).split('.fasta')[0] #make the bin_name the filename minus the word fasta
    for record in SeqIO.parse(f, 'fasta'): #using SeqIO package look for the item name then do
        results.append((record.id, bin_name))#add the results of the parse and only take the record.id (contig name) and the bin_name from above

process_file(results, 'maxbin.bins.tsv') #make the file from all the results from the results.append function into <filename>
#metabat

results = [] #clear the results object 
file_list = glob.glob('./metabat2/*.fa')
for f in file_list:
    bin_name = os.path.basename(f).split('.fa')[0]
    for record in SeqIO.parse(f, 'fasta'):
        results.append((record.id, bin_name))


process_file(results, 'metabat.bins.tsv')        

#binsanity
results = []
file_list = glob.glob('./binsanity/IGM-BinsanityWF/*.fna')
for f in file_list:
    bin_name = os.path.basename(f).split('.fna')[0]
    bin_name = bin_name.replace('cellular.contigs.from.anvio_', '')
    for record in SeqIO.parse(f, 'fasta'):
        results.append((record.id, bin_name))

process_file(results, 'binsanity.bins.tsv')
```
#anvi'o import a collection and export it allows for it to be exported as misc data to allow for viewing my multiple binning softwares
#metabat2
anvi-import-collection \
metabat2.bins.tsv \
-p SAMPLES-MERGED/PROFILE.db \
-c contigs.db \
-C MetaBat2 \
--contigs-mode

anvi-export-collection \
-p SAMPLES-MERGED/PROFILE.db \
-C MetaBat2 \
-O MetaBat2

anvi-import-misc-data \
MetaBat2.txt \
-p SAMPLES-MERGED/PROFILE.db \
-t items

#binsanity
anvi-import-collection \
binsanity.bins.tsv \
-p SAMPLES-MERGED/PROFILE.db \
-c contigs.db \
-C BinSanity \
--contigs-mode

anvi-export-collection \
-p SAMPLES-MERGED/PROFILE.db \
-C BinSanity \
-O BinSanity

anvi-import-misc-data \
BinSanity.txt \
-p SAMPLES-MERGED/PROFILE.db \
-t items

#MaxBin
anvi-import-collection \
maxbin.bins.tsv \
-p SAMPLES-MERGED/PROFILE.db \
-c contigs.db \
-C MaxBin \
--contigs-mode

anvi-export-collection \
-p SAMPLES-MERGED/PROFILE.db \
-C MaxBin \
-O MaxBin

anvi-import-misc-data \
MaxBin.txt \
-p SAMPLES-MERGED/PROFILE.db \
-t items

#run anvi'o plot
anvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db
