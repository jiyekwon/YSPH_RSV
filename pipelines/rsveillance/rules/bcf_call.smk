
rule bcf_call:
    input:
        tocall = lambda wildcards: expand('results/ivar/{sample}-{{target}}-itrim.bam',sample=get_mash_samples(wildcards)),
        indexed = lambda wildcards: expand("results/ivar/{sample}-{{target}}-itrim.bam.bai",sample=get_mash_samples(wildcards)),
    output:
        vcf=temporary('results/bcftools/{target}-all-unfilt.vcf.gz')
    params:
        ref=os.path.join(config['refsdir'],'{target}.fasta')
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        stderr="logs/ivar/bcf_call-{target}.err"
    message: "Calling variants for all samples"
    shell:
        """
        bcftools mpileup -Ou -o variants.bcf -f {params.ref} {input.tocall} 2>&1  > {log.stderr}
        bcftools call --ploidy 1 -vcO z -o {output.vcf} variants.bcf 2>&1  >> {log.stderr}
        tabix -p vcf {output.vcf} 2>&1  >> {log.stderr}	        
	"""

rule bcf_filter:
    input:
        vcf='results/bcftools/{target}-all-unfilt.vcf.gz'
    output:
        vcf_filt='results/bcftools/{target}-all.vcf.gz',
        vcf_filt_idx='results/bcftools/{target}-all.vcf.gz.tbi'
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

