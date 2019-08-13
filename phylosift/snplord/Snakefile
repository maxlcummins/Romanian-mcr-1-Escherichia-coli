configfile:
	"config_files/snake-snippy.config.json"

sample_ids, = glob_wildcards(config['reads']+"/{sample}.R1.fastq.gz")
REF, = glob_wildcards(config['ref']+"/{reference}.fa")


print(sample_ids,)
print(REF,)

rule all:
	input:
		expand("data/output/snippyout/{reference}.{sample}.out", sample=sample_ids, reference=REF),		
		expand("data/output/fasttree/{reference}.clean.fullcore.tree", reference=REF),
		expand("data/output/snp_dists/{reference}.pairwise_snps.csv", reference=REF)		

rule snippy_run:
	input:
		ref = config['ref']+"/{reference}.fa",
		r1 = config['reads']+"/{sample}.R1.fastq.gz",
		r2 = config['reads']+"/{sample}.R2.fastq.gz"
	output: 
		directory("data/output/snippyout/{reference}.{sample}.out")
	conda: 
		"config_files/snippy.yaml"
	shell: 
		"snippy --outdir {output} --ref {input.ref} --R1 {input.r1} --R2 {input.r2}"


rule snippy_core:
	input:
		ref = config['ref']+"/{reference}.fa",	
	output:	
		"data/output/core/{reference}.full.aln"
	conda:
		"config_files/snippy.yaml"
	shell:
		"snippy-core --prefix data/output/core/{wildcards.reference} --ref {input.ref} data/output/snippyout/*.out"

rule snippy_clean:
	input:
		"data/output/core/{reference}.full.aln"
	output:
		"data/output/core/{reference}.clean.full.aln"
	conda: 
		"config_files/snippy.yaml"	
	shell:
		"snippy-clean_full_aln {input} > {output}"

rule gubbarnt:
	input:
		"data/output/core/{reference}.clean.full.aln"
	output:
		"data/output/gubbins/gubbins.{reference}.filtered_polymorphic_sites.fasta"
	params:
		prefix = "data/output/gubbins/gubbins.{reference}"
	conda:
		"config_files/gubbins.yaml"
	shell:
		"""
		run_gubbins.py -v -p {params.prefix} {input}
		rm {wildcards.reference}.clean.full.aln.seq.joint.txt
		"""
rule snp_sites:
	input:
		"data/output/gubbins/gubbins.{reference}.filtered_polymorphic_sites.fasta"
	output:
		"data/output/snp_sites/{reference}.clean.fullcore.aln"
	conda: 
		"config_files/snippy.yaml"
	shell:
		"snp-sites -c {input} > {output}"

rule snp_dists:
	input:
		"data/output/snp_sites/{reference}.clean.fullcore.aln"
	output:
		"data/output/snp_dists/{reference}.pairwise_snps.csv"
	conda: 
		"config_files/snp_dists.yaml"
	shell: 
		"snp-dists -c {input} > {output}"

rule fasttree:
	input:
		"data/output/snp_sites/{reference}.clean.fullcore.aln"
	output:
		"data/output/fasttree/{reference}.clean.fullcore.tree"
	conda:
		"config_files/fasttree.yaml"
	shell:
		"fasttree -gtr -nt {input} > {output}"






