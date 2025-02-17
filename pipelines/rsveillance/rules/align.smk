
rule bwa_align_pipeline:
    input:
        bam = 'results/align/{sample}_{target}.bam',
        stat = 'results/align/{sample}_{target}.flagstat',
        depth = 'results/align/{sample}_{target}.depth',
        subsamp = 'results/align/{sample}_{target}.substat'


rule bwa_index:
    input:
        fasta = "results/refs/{target}.fasta"
    output:
        bwt = "results/refs/{target}.fasta.bwt"
    log:
        log = "logs/bwa/index_{target}.log"
    params:
        idx_script = "scripts/indexer.sh",
    resources:
        mem_mb="40G",
        cpus_per_task=4,
        runtime=300
    container: "docker://sethnr/pgcoe_anypipe:0.01"
    shell:"""
        {params.idx_script} \
                -C {resources.cpus_per_task} \
                {input.fasta} >> {log.log} 2>&1
            """


rule bwa_align:
    input:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz',
        indexedref="results/refs/{target}.fasta.bwt"
    output:
        aligned = temporary('results/align/{sample}_{target}_unsort.bam'),
        #aligned = 'results/align/{sample}_{target}_unsort.bam',
    group: "aligngroup"
    params:
        ref="results/refs/{target}.fasta",
    resources:
        runtime=180,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 8000),
        cpus_per_task=4,
        cores=1,
    log:
        stderr="logs/align/bwa_mem_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        echo "Aligning reads for {wildcards.sample} to {params.ref}\n"
        echo 'bwa mem -t {resources.cpus_per_task} {params.ref} {input.R1} {input.R2} \n  
                 -o {output.aligned}  \n' 
        bwa mem -t {resources.cpus_per_task} {params.ref} {input.R1} {input.R2} \
            -o {output.aligned}  2> {log.stderr}
        """

rule flagstat:
    input:
        aligned = 'results/align/{sample}_{target}_unsort.bam',
    output:
        flagstat = 'results/align/{sample}_{target}.flagstats'
    params:
        ref="results/refs/{target}.fasta",
        sleeplen=60,
    group: "aligngroup"
    resources:
        runtime=30,
        mem_mb=4000,
        cores=1,
    log:
        stderr="logs/align/flagstat_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        samtools quickcheck {input.aligned} 2>&1 > {log.stderr}
        if [[ $? >0 ]]; then
            echo "sleeping {params.sleeplen}"
            sleep {params.sleeplen}
        else:
           echo "quickcheck good $?"
        fi

        samtools flagstat -@ {resources.cores} -O tsv   {input.aligned} > {output.flagstat}

        """

rule remove_unaligned:
    input:
        unsorted = 'results/align/{sample}_{target}_unsort.bam',
    output:
        aligned = temporary('results/align/{sample}_{target}_aligned.bam'),
    params:
        ref="results/refs/{target}.fasta",
        sleeplen=60,
    group: "aligngroup"
    resources:
        runtime=30,
        mem_mb=4000,
        cores=1,
    log:
        stderr="logs/align/flagstat_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        samtools view -b -F 4 -F 2048 -o {output.aligned} {input.unsorted}
        """

rule sam_subsample:
    input:
        aligned = 'results/align/{sample}_{target}_aligned.bam'
    output:
        subsamp = temporary('results/align/{sample}_{target}_thin.bam'),
        subfactor = 'results/align/{sample}_{target}.substat'
    resources:
        mem_mb=8000,
        runtime=60,
        cores=4,
    group: "aligngroup"
    params:
        maxsize_mb=1024
    log:
        stderr="logs/align/{sample}_{target}_sort.err"
    run:
        import subprocess
        insize = round(os.stat(input.aligned).st_size / 1048576)
        subfact = get_subsample_factor(insize, params.maxsize_mb)
        subfraction = 1/subfact
        if subfact < 1:
            subprocess.run(["samtools", "view","-b","-s",str(subfraction),"-o",output.subsamp,input.aligned])
        else:
            subprocess.run(["cp",input.aligned,output.subsamp])
        #subprocess.run(["samtools","index",output.subsamp])
        f = open(output.subfactor, "w")
        print("\t".join([wildcards.sample,str(subfact)]),file=f)
        f.close()

rule sam_sort:
    input:
        aligned = 'results/align/{sample}_{target}_aligned.bam'
    output:
        sorted = 'results/align/{sample}_{target}.bam',
        idx = 'results/align/{sample}_{target}.bam.bai'
    resources:
        mem_mb=16000,
        runtime=40,
        cores=4,
    group: "aligngroup"
    log:
        stderr="logs/align/sort_{sample}_{target}.err"
    message: "Sorting and indexing reads"
    shell:
        """
        samtools sort -@ {resources.cores} -o {output.sorted} {input.aligned} 2>&1 > {log.stderr}
        samtools index {output.sorted} 2>&1 >> {log.stderr}
        """



