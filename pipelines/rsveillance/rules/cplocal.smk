rule cp_local_fq:
    input:
        read_location=os.path.join(config['readdir'],'{sample}'),
    output:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz'
    log:
        log = "logs/datacopy/{sample}.log"
    group: "cplocal"
    resources:
        mem_mb="2G",
        cpus_per_task=1,
        runtime=30
    container: None
    shell:"""
        #find {input.read_location}

        R1=` find {input.read_location}/ | grep _R1_ `
        R2=` find {input.read_location}/ | grep _R2_ `

        echo cp $R1 {output.R1}
        cp $R1 {output.R1}
        
        echo cp $R2 {output.R2}
        cp $R2 {output.R2}
        
        """


rule cp_local_data:
    input:
        refdir=config['refsdir'],
        ampdir=config['ampsdir'],
        refsfile=config['refs']
    output:
        fa = expand('results/refs/{target}.fasta',target=TARGETS)
    log:
        log = "logs/datacopy/copyrefs.log"
    resources:
        mem_mb="2G",
        cpus_per_task=1,
        runtime=60
    params:
        targets = TARGETS
    container: None
    shell:"""
        cp -LR {input.refdir}/* results/refs/
        cp -LR {input.ampdir}/* results/refs/
        cp {input.refsfile} results/refs/refs.txt        
        """

