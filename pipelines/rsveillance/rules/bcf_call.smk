
rule bcf_call:
    input:
        tocall = lambda wildcards: expand("results/bams/{{sample}}_{target}_itrim.bam",target=get_mash_calls(wildcards)),
        indexed = lambda wildcards: expand("results/bams/{{sample}}_{target}_itrim.bam.bai",target=get_mash_calls(wildcards)),
        #tocall=expand('results/ivar/{sample}_{target}_itrim.bam',sample=SAMPLES),
        #indexed=expand('results/ivar/{sample}_{target}_itrim.bam.bai'),sample=SAMPLES)
    output:
        vcf=temporary('results/bcftools/{target}_all_unfilt.vcf.gz')
    params:
        ref='{refs}/{target}.fa'
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        stderr="logs/ivar/bcf_call.err"
    message: "Calling variants for all samples"
    shell:
        """
        bcftools mpileup -Ou -o variants.bcf -f {params.ref} {input.tocall} 1> {log.stdout} 2> {log.stderr}
        bcftools call --ploidy 1 -vcO z -o {output.vcf} variants.bcf 1>> {log.stdout} 2>> {log.stderr}
        tabix -p vcf {output.vcf} 1>> {log.stdout} 2>> {log.stderr}	        
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
        stderr="logs/bcf_filter.err"
    shell:
        """
        bcftools filter -O z -o {output.vcf_filt} -i 'QUAL>10 & DP>10' {input.vcf} 2>&1 > {log.stderr}
        """

