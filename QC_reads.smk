configfile: "configs/sample_configs.yaml"

rule all:
	input:
		expand("{sample}.qtrimmed.fq.gz", sample=config["sample"]),
		expand("{sample}.merged.fq.gz", sample=config["sample"]),
		directory(expand("{sample}_spades_out/", sample=config["sample"]))

rule rm_duplicates:
	input: "{sample}.fq.gz"
	output: temp("{sample}.clumped.fq.gz")
	shell: "clumpify.sh in={input} out={output} dedupe optical"

rule rm_low_qual_regions:
	input: temp("{sample}.clumped.fq.gz")
	output: temp("{sample}.filtered_by_tile.fq.gz")
	shell: "filterbytile.sh in={input} out={output}"

rule trim_adaptors:
	input: temp("{sample}.filtered_by_tile.fq.gz")
	output: temp("{sample}.trimmed.fq.gz")
	shell: "bbduk.sh in={input} out={output} ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered"

rule rm_synth_artifacts:
	input: temp("{sample}.trimmed.fq.gz")
	output: temp("{sample}.filtered.fq.gz")
	shell: "bbduk.sh in={input} out={output} k=31 ref=artifacts,phix ordered cardinality"

rule err_cor_1:
	input: temp("{sample}.filtered.fq.gz")
	output: temp("{sample}.ecco.fq.gz")	
	shell: "bbmerge.sh in={input} out={output} ecco mix vstrict ordered"

rule err_cor_2:
	input: temp("{sample}.ecco.fq.gz")
	output: temp("{sample}.eccc.fq.gz")
	shell: "clumpify.sh in={input} out={output} ecc passes=4 reorder"

rule err_cor_3:
	input: temp("{sample}.eccc.fq.gz")
	output: temp("{sample}.ecct.fq.gz")
	shell: "tadpole.sh in={input} out={output} ecc k=62 ordered prefilter=2"

rule normalise:
	input: temp("{sample}.ecct.fq.gz")
	output: temp("{sample}.normalized.fq.gz")
	shell: "bbnorm.sh in={input} out={output} target=100"

rule overlapping_reads:
	input: temp("{sample}.normalized.fq.gz")
	output: 
		merged="{sample}.merged.fq.gz",
		unmerged=("{sample}.unmerged.fq.gz")
	shell: "bbmerge-auto.sh in={input} out={output.merged} outu={output.unmerged} strict k=93 extend2=80 rem ordered prefilter=2"

rule QC_unmerged:
	input: temp("{sample}.unmerged.fq.gz")
	output: "{sample}.qtrimmed.fq.gz"
	shell: "bbduk.sh in={input} out={output} qtrim=r trimq=10 minlen=70 ordered"

rule SPAdes_assembly:
	input: 
		merged="{sample}.merged.fq.gz",
		qtrimmed="{sample}.qtrimmed.fq.gz"
	params: meta="--meta"
	output: directory("{sample}_spades_out")
	shell: "spades.py {params.meta} -k 25,55,95,125 --phred-offset 33 -s {input.merged} --12 {input.qtrimmed} -o {output}"

