
rule bcf_call:
    input:
        tocall = lambda wildcards: expand('results/ivar/{sample}_{{target}}_itrim.bam',sample=get_mash_samples(wildcards)),
        indexed = lambda wildcards: expand("results/ivar/{sample}_{{target}}_itrim.bam.bai",sample=get_mash_samples(wildcards)),
        ref="results/refs/{target}.fasta"
    output:
        bamlist=temporary('results/bcftools/bamlist_{target}.txt'),
        bcf=temporary('results/bcftools/{target}_variants.vcf'),
        vcf=temporary('results/bcftools/{target}_all_unfilt.vcf.gz')
    resources:
        mem_mb=24000,
        runtime=240,
        cpus_per_task=4,
    params:
        threads=3
    log:
        stderr="logs/ivar/bcf_call_{target}.err"
    message: "Calling variants for all samples"
    shell:
        """
        ls -1 results/ivar/*_{wildcards.target}_itrim.bam > {output.bamlist}
        bcftools mpileup -Ov --threads {params.threads} -o {output.bcf} -f {input.ref} -b {output.bamlist} #2>&1  > {log.stderr}
        bcftools call --threads {params.threads} --ploidy 1 -A -vcO z -o {output.vcf} {output.bcf}         #2>&1  >> {log.stderr}
        tabix -p vcf {output.vcf} 2>&1  >> {log.stderr}	        
	"""

rule bcf_filter:
    input:
        vcf='results/bcftools/{target}_all_unfilt.vcf.gz'
    output:
        vcf_filt='results/bcftools/{target}_all.vcf.gz',
        vcf_filt_idx='results/bcftools/{target}_all.vcf.gz.tbi'
    resources:
        mem_mb=8000,
        runtime=180,
    log:
        stderr="logs/bcf_filt_{target}.err"
    shell:
        """
        bcftools filter -O z -o {output.vcf_filt} -i 'QUAL>10 & DP>10' {input.vcf} 2>&1 > {log.stderr}
	    tabix -p vcf {output.vcf_filt}
        """



rule bcf_call_untrim:
    input:
        tocall = lambda wildcards: expand('results/align/{sample}_{{target}}.bam',sample=get_mash_samples(wildcards)),
        indexed = lambda wildcards: expand("results/align/{sample}_{{target}}.bam.bai",sample=get_mash_samples(wildcards)),
        ref="results/refs/{target}.fasta"
    output:
        bamlist=temporary('results/bcftools/bamlist_untrim_{target}.txt'),
        bcf=temporary('results/bcftools/{target}_untrim_variants.vcf'),
        vcf=temporary('results/bcftools/{target}_untrim_all_unfilt.vcf.gz')
    resources:
        mem_mb=16000,
        runtime=240,
        cpus_per_task=4,
    params:
        threads=3
    log:
        stderr="logs/ivar/bcf_call_untrim_{target}.err"
    message: "Calling variants for all samples"
    shell:
        """
        ls -1 results/align/*_{wildcards.target}.bam > {output.bamlist}
        bcftools mpileup -Ov --threads {params.threads} -o {output.bcf} -f {input.ref} -b {output.bamlist} 2>&1  > {log.stderr}
        bcftools call --threads {params.threads} --ploidy 1 -A -vcO z -o {output.vcf} {output.bcf} 2>&1  >> {log.stderr}
        tabix -p vcf {output.vcf} 2>&1  >> {log.stderr}	        
	"""

rule bcf_filter_untrim:
    input:
        vcf='results/bcftools/{target}_untrim_all_unfilt.vcf.gz'
    output:
        vcf_filt='results/bcftools/{target}_untrim_all.vcf.gz',
        vcf_filt_idx='results/bcftools/{target}_untrim_all.vcf.gz.tbi'
    resources:
        mem_mb=8000,
        runtime=180,
    log:
        stderr="logs/bcf_filt_{target}.err"
    shell:
        """
        bcftools filter -O z -o {output.vcf_filt} -i 'QUAL>10 & DP>10' {input.vcf} 2>&1 > {log.stderr}
	    tabix -p vcf {output.vcf_filt}
        """

