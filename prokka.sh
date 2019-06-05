prokka $i \
--outdir $(basename $i .reformatted.fasta) \
--force \
--prefix $(basename $i .reformatted.fasta) \
--addgenes \
--addmrna \
--locustag gene \
--compliant \
--genus pelagibacter \
--species ubique \
--strain $(basename $i .reformatted.fasta) \
--cpus 0 \
--rfam
