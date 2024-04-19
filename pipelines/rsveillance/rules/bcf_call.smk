
rule bcf_call:
    input:
        tocall = lambda wildcards: expand('results/bams/{{sample}}_{target}_itrim.bam',target=get_mash_targets(wildcards)),
        indexed = lambda wildcards: expand("results/bams/{{sample}}_{target}_itrim.bam.bai",target=get_mash_targets(wildcards)),
        #tocall=expand('results/ivar/{sample}_{target}_itrim.bam',sample=SAMPLES),
    output:
        vcf=temporary('results/bcftools/{target}_all_unfilt.vcf.gz')
    params:
        ref=config['refsdir']+'/{target}.fa'
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        stderr="logs/ivar/bcf_call_{target}.err"
    message: "Calling variants for all samples"
    shell:
        """
        bcftools mpileup -Ou -o variants.bcf -f {params.ref} {input.tocall} 2>&1  > {log.stderr}
        bcftools call --ploidy 1 -vcO z -o {output.vcf} variants.bcf 2>&1  >> {log.stderr}
        tabix -p vcf {output.vcf} 2>&1  >> {log.stderr}	        
	"""

rule bcf_filter:
    input:
        vcf='results/bcftools/{target}_all_unfilt.vcf.gz'
    output:
        vcf_filt='results/bcftools/{target}_all_filt.vcf.gz'
    resources:
        mem_mb=8000,
        runtime=180,
    log:
        stderr="logs/bcf_filt_{target}.err"
    shell:
        """
        bcftools filter -O z -o {output.vcf_filt} -i 'QUAL>10 & DP>10' {input.vcf} 2>&1 > {log.stderr}
        """

