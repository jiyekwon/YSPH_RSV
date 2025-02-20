
rule bcf_call_one:
    input:
        tocall = 'results/ivar/{sample}_{target}_itrim.bam',
        indexed = "results/ivar/{sample}_{target}_itrim.bam.bai",
        ref="results/refs/{target}.fasta"
    output:
        pileup=temporary('results/recondata/{sample}_{target}_pileup.vcf'),
        unfilt=temporary('results/recondata/{sample}_{target}_unfilt.vcf')
        vcf=temporary('results/recondata/{target}/vcf/{sample}.vcf')
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
        bcftools mpileup -Ov --threads {params.threads} -o {output.bcf} -f {input.ref} {input.tocall}
        bcftools call --threads {params.threads} --ploidy 1 -A -vcO v -o {output.unfilt} {output.pileup}
        bcftools filter -o {output.vcf} -i 'QUAL>10 & DP>10' {output.unfilt}
	"""



rule align_fastas:
    input:
        ref="results/refs/{target}.fasta",
        callfile='results/summary/final_calls.txt',
        fastas = lambda wildcards: expand('results/ivar/{sample}_consensus.fa',sample=get_called_target_samples(wildcards)),
    output:
        consensus='results/recondata/{target}_consensus.fasta',
        aligned = 'results/recondata/{target}/aligned.fasta',
    shell:
    """
        cat {input.fastas} > {output.consensus}

        mafft --6merpair --keeplength --addfragments {input.consensus} {target.ref} > {output.aligned}
    """



rule collect_files:
    input:
        ref = 'results/recondata/{target}/aligned.fasta',
        aligned = 'results/recondata/{target}/aligned.fasta',
        depth = 'results/recondata/{target}/depth.txt',
        date
        vcfs = lambda wildcards: expand('results/recondata/{target}/vcfs/{sample}.vcf',sample=get_called_target_samples(wildcards,input.callfile)),
        callfile='results/summary/final_calls.txt',        
    output:
        finished=temporary('results/recondata/{target}/collected.txt')
    shell:
        """ 
        touch {output.finished}
        """






def get_called_target_samples(wildcards,callfile):
    print("mash: requesting samples for "+wildcards.target+" from "+callfile,file=sys.stderr)  
    with callfile.open() as f:
        mashlines = f.read().splitlines()
    mysamples = [L.split()[0] for L in mashlines if L.split()[2]==wildcards.target]
    return mysamples
