
rule ivar_pipeline:
    input:
        consensus='results/ivar/{sample}_{target}_consensus.fa',
        ivariants='results/ivar/{sample}_{target}_ivariants.tsv'


rule ivar_pclip:
    input:
        aligned = 'results/align/{sample}_{target}.bam',
        primers = '{refs}/{target}.bed'
    output:
        aln_trimmed='results/ivar/{sample}_{target}_itrim.bam'
    resources:
            mem_mb=8000,
            runtime=1440,
    log:
        stdout="logs/ivar/{sample}_trim.out",
        stderr="logs/ivar/{sample}_trim.err"
    message: "QC and soft-clipping primers using iVar"
    shell:
        """
        ivar trim -i {input.aligned} -b {input.primers} -p {output.aln_trimmed} -e 2>&1 > {log.stderr}
        """


rule sam_pileup:
    input:
        tocall='results/ivar/{sample}_{target}_itrim.bam',
        indexed='results/ivar/{sample}_{target}_itrim.bam.bai'
    output:
        pileup=temporary('results/ivar/{sample}_{target}.pileup'),
    params:
        ref=config['reference'],
        threshold=0.2,
        depth=20,
    log:
        stderr="logs/ivar/pileup_{sample}_{target}.err"
    shell:
        """
        samtools mpileup -aa -A -d 0 -Q 0 -f {params.ref} -o {output.pileup} {input.tocall} 2>&1 > {log.stderr}
        """

rule ivar_consensus:
    input:
        pileup='results/ivar/{sample}_{target}.pileup'
    output:
        consensus='results/ivar/{sample}_{target}_consensus.fa'
    params:
        ref=config['reference'],
        threshold=0.2,
        depth=20,
        prefix='results/ivar/{sample}_consensus'
    log:
        stderr="logs/ivar/consensus_{sample}_{target}.err",
    shell:
        """
        cat {input.pileup} | ivar consensus -t {params.threshold} \
           -m {params.depth} -p {params.prefix} -i {wildcards.sample} \
           2>&1 > {log.stderr}
        """

rule ivar_variants:
    input:
        pileup='results/ivar/{sample}_{target}.pileup'
    output:
        ivariants='results/ivar/{sample}_{target}_ivariants.tsv',
    params:
        ref=config['reference'],
        threshold=0.2,
        depth=20,
        qual=20,
        prefix='results/ivar/{sample}_{target}_ivariants'
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        stderr="logs/ivar/{sample}_consensus.err"
    shell:
        """
        cat {input.pileup} | \
         ivar variants -q {params.qual} -r {params.ref} -t {params.threshold} \
           -m {params.depth} -p {params.prefix}  2>&1 >  {log.stderr}        
        """
