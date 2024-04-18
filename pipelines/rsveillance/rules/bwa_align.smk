
rule bwa_align_pipeline:
    input:
        bam = 'results/align/{sample}_{target}.bam',
        stat = 'results/align/{sample}_{target}.flagstat',
        depth = 'results/align/{sample}_{target}.depth',
        subsamp = 'results/align/{sample}_{target}.substat'

rule bwa_index:
    input:
        fasta = config['refsdir']+"/{target}.fasta"
    output:
        bwt = config['refsdir']+"/{target}.fasta.bwt"
    log:
        log = "/logs/bwa/index_{target}.log"
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
        read_location=config['readdir']+'/{sample}',
        indexedref=config['readdir']+"/{target}.fasta.amb"
    params:
        ref=config['refsdir']+"/{target}.fasta"
    resources:
        mem_mb=8000,
        runtime=240,
    output:
        aligned = temporary('results/align/{sample}_{target}_unsort.bam')
    log:
        stderr="logs/align/bwa_mem_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        echo "Aligning reads for {wildcards.sample} to {params.ref}\n"
        echo 'bwa mem -o {output.aligned} {params.ref} {input.read_location}/*R1* {input.read_location}/*R2* \n' 
        bwa mem  -o {output.aligned} {params.ref} {input.read_location}/*R1* {input.read_location}/*R2*  >> {log.stdout} 2>> {log.stderr}
        """

rule sam_sort:
    input:
        aligned = 'results/align/{sample}_{target}_unsort.bam'
    output:
        aligned = temporary('results/align/{sample}_{target}_sort.bam'),
        idx = 'results/align/{sample}_{target}_sort.bam.idx'
    resources:
        mem_mb=16000,
        runtime=180,
	cores=4,
    log:
        stderr="logs/align/sort_{sample}_{target}.err"
    message: "Sorting and indexing reads"
    shell:
        """
        samtools view -b -F 4 -F 2048 {input.aligned} | samtools sort -@ {resources.cores} -o {output.aln_itrim_sorted} 2>&1 > {log.stderr}
        samtools index {output.aligned} 2>&1 >> {log.stderr}
        """

rule sam_subsample:
    input:
        aligned = 'results/align/{sample}_{target}_sort.bam'
    output:
        subsamp = 'results/align/{sample}_{target}.bam',
        idx = 'results/align/{sample}_{target}.bam.idx',
        subfactor = temporary('results/align/{sample}_{target}.substat')
    resources:
        mem_mb=16000,
        runtime=180,
    params:
        subfactor=lambda getfact: get_subsample_factor(input.size_mb,1024) , 
        cores=4,
    log:
        stderr="logs/align/{sample}_{target}_sort.err"
    shell:
        """
        samtools view -b -s {params.subfactor} {input.aln_itrim_sorted} -o {output.subsamp} 2>&1 > {log.stderr}
        samtools index {output.subsamp}

        echo "{params.subfactor} subfactor\n" > {output.subfactor}
        """


rule indexbam:
    input:
        bam = '{samplename}.bam'
    output:
        indexedbam = '{samplename}.bam.bai'
    log:
        stdout="logs/bamindex/{samplename}.out",
        stderr="logs/bamindex/{samplename}.err"
    shell:
        """
        samtools index {input.bam} > {log.stdout} 2> {log.stderr}
        """




### STATS ###

rule depth:
    input:
        bams='results/align/{sample}_{target}.bam'
    output:
        'results/align/{sample}_{target}_depth.txt'
    resources:
        mem_mb=8000,
        runtime=180,
    params:
        maxdepth=0,
        minmapqual=60,
        minbasequal=13
    log:
        stderr="logs/depth/{sample}_{target}.err"
    shell:
        """
        samtools depth -a -H {input.bams} -o {output} 2>&1 >  {log.stderr}
        """

rule flagstat:
    input:
        bam='results/align/{sample}_{target}.bam'
    output:
        'results/align/{sample}_{target}_flagstat.txt'
    resources:
        mem_mb=8000,
        runtime=180,
    params:
        maxdepth=0,
        minmapqual=60,
        minbasequal=13
    log:
        stderr="logs/align/{sample}_{target}_flagstat.err"
    shell:
        """
        samtools flagstats -O tsv {input.bam} 1> {output} 2> {log.stderr}
        """

