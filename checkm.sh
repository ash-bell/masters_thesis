checkm taxon_set class Alphaproteobacteria Alphaproteobacteria_markers
checkm analyze Alphaproteobacteria_markers <folder_with_SAGs_assemblies> checkm_output -x fasta -t 16
checkm qa --out_format 2 --file checkm_output/completeness_stats.txt --tab_table -t 16 Alphaproteobacteria_markers checkm_output
