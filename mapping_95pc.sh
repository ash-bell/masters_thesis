bowtie2-build <genome> <genome_name>
bowtie2 \
--threads 16 \
-x <genome_name> \
-1 <fwd_read> \
-2 <rev_read> \
--no-unal \
-S <output_sam_file>

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
