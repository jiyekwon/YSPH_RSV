
rule bwa_align_pipeline:
    input:
        bam = 'results/align/{sample}-{target}.bam',
        stat = 'results/align/{sample}-{target}.flagstat',
        depth = 'results/align/{sample}-{target}.depth',
        subsamp = 'results/align/{sample}-{target}.substat'


rule bwa_index:
    input:
        fasta = os.path.join(config['refsdir'],"{target}.fasta")
    output:
        bwt = os.path.join(config['refsdir'],"{target}.fasta.bwt")
    log:
        log = "logs/bwa/index-{target}.log"
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
        read_location=os.path.join(config['readdir'],'{sample}'),
        indexedref=os.path.join(config['refsdir'],"{target}.fasta.bwt")
    output:
        aligned = temporary('results/align/{sample}-{target}-unsort.bam')
    params:
        ref=config['refsdir']+"{target}.fasta"
    resources:
        runtime= 600,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 4000),
	cores=4,
        #mem_mb=8000,
	#runtime=lambda _: max(60* (input.size_mb/1024), 120),
    log:
        stderr="logs/align/bwa_mem_{sample}-{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        echo "Aligning reads for {wildcards.sample} to {params.ref}\n"
        echo 'bwa mem -o {output.aligned} {params.ref} {input.read_location}/*R1* {input.read_location}/*R2* \n' 
        bwa mem -t {resources.cores} {params.ref} {input.read_location}/*R1* {input.read_location}/*R2* | \
            samtools view -b -F 4 -F 2048 1> {output.aligned}  2> {log.stderr}
        """

rule sam_subsample:
    input:
        aligned = 'results/align/{sample}-{target}-unsort.bam'
    output:
        subsamp = temporary('results/align/{sample}-{target}-thin.bam'),
        #idx = 'results/align/{sample}-{target}-thin.bam.bai',
        #subfactor = temporary('results/align/{sample}-{target}.substat')	
        subfactor = 'results/align/{sample}-{target}.substat'
    resources:
        mem_mb=16000,
        runtime=180,
        cores=4,
    params:
        subfactor=lambda _: get_subsample_factor( 2048 ) ,
        #subfactor=lambda _: get_subsample_factor( input.size_mb, 1024) ,
        #subfactor=0.5
	maxsize_mb=1024
    log:
        stderr="logs/align/{sample}-{target}-sort.err"
#    shell:
#        """
#        samtools view -b -s {params.subfactor} {input.aligned} -o {output.subsamp} 2>&1 > {log.stderr}
#        samtools index {output.subsamp}
#
#        echo "{params.subfactor} subfactor\n" > {output.subfactor}
#        """
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
        print(subfact,file=f)
        f.close()

rule sam_sort:
    input:
        aligned = 'results/align/{sample}-{target}-thin.bam'
    output:
        sorted = 'results/align/{sample}-{target}.bam',
        idx = 'results/align/{sample}-{target}.bam.bai'
    resources:
        mem_mb=16000,
        runtime=180,
	cores=4,
    log:
        stderr="logs/align/sort-{sample}-{target}.err"
    message: "Sorting and indexing reads"
    shell:
        """
        samtools sort -@ {resources.cores} -o {output.sorted} {input.aligned} 2>&1 > {log.stderr}
        samtools index {output.sorted} 2>&1 >> {log.stderr}
        """




rule indexbam:
    input:
        bam = '{samplename}.bam'
    output:
        indexedbam = '{samplename}.bam.bai'
    log:
        stderr="logs/bamindex/{samplename}.err"
    shell:
        """
        samtools index {input.bam} 2>&1 >  {log.stderr}
        """




### STATS ###

rule depth:
    input:
        bams='results/align/{sample}-{target}.bam'
    output:
        depth='results/align/{sample}-{target}-depth.txt'
    resources:
        mem_mb=8000,
        runtime=180,
    params:
        maxdepth=0,
        minmapqual=60,
        minbasequal=13
    log:
        stderr="logs/depth/{sample}-{target}.err"
    shell:
        """
        samtools depth -a -H {input.bams} -o {output.depth} 2>&1 >  {log.stderr}
        """

rule flagstat:
    input:
        bam='results/align/{sample}-{target}.bam'
    output:
        'results/align/{sample}-{target}-flagstat.txt'
    resources:
        mem_mb=8000,
        runtime=180,
    params:
        maxdepth=0,
        minmapqual=60,
        minbasequal=13
    log:
        stderr="logs/align/{sample}-{target}-flagstat.err"
    shell:
        """
        samtools flagstats -O tsv {input.bam} 1> {output} 2> {log.stderr}
        """

