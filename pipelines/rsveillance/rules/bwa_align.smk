
rule bwa_align_pipeline:
    input:
        bam = 'results/align/{sample}_{target}.bam',
        stat = 'results/align/{sample}_{target}.flagstat',
        depth = 'results/align/{sample}_{target}.depth',
        subsamp = 'results/align/{sample}_{target}.substat'


rule bwa_index:
    input:
        fasta = os.path.join(config['refsdir'],"{target}.fasta")
    output:
        bwt = os.path.join(config['refsdir'],"{target}.fasta.bwt")
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
        runtime=60
    container: None
    shell:"""
        find {input.read_location}

        R1=` find {input.read_location} | grep _R1_ `
        R2=` find {input.read_location} | grep _R2_ `

        echo cp $R1 {output.R1}
        cp $R1 {output.R1}
        
        echo cp $R2 {output.R2}
        cp $R2 {output.R2}
        
        """


rule bwa_align:
    input:
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz',
        indexedref=os.path.join(config['refsdir'],"{target}.fasta.bwt")
    output:
        #aligned = temporary('results/align/{sample}_{target}_unsort.bam'),
        aligned = 'results/align/{sample}_{target}_unsort.bam',
    group: "aligngroup"
    params:
        ref=config['refsdir']+"{target}.fasta"
    resources:
        runtime= 600,
        mem_mb=lambda wc, input: max(2 * input.size_mb, 4000),
	    cores=4,
    log:
        stderr="logs/align/bwa_mem_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        echo "Aligning reads for {wildcards.sample} to {params.ref}\n"
        echo 'bwa mem -o {output.aligned} {params.ref} {input.R1} {input.R2} \n' 
        bwa mem -t {resources.cores} {params.ref} {input.R1} {input.R2} \
            1> {output.aligned}  2> {log.stderr}
        """

rule flagstat:
    input:
        aligned = 'results/align/{sample}_{target}_unsort.bam',
    output:
        flagstat = 'results/align/{sample}_{target}.flagstats'
    params:
        ref=config['refsdir']+"{target}.fasta",
        sleeplen=60,
    group: "aligngroup"
    resources:
        runtime= 600,
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

rule flagstat:
    input:
        unsorted = 'results/align/{sample}_{target}_unsort.bam',
    output:
        aligned = 'results/align/{sample}_{target}_aligned.bam',
    params:
        ref=os.path.join(config['refsdir'],"{target}.fasta"),
        sleeplen=60,
    group: "aligngroup"
    resources:
        runtime= 600,
        mem_mb=4000,
	    cores=1,
    container: "docker://sethnr/aaaa_gvs_processing:0.02"
    log:
        stderr="logs/align/flagstat_{sample}_{target}.err",
    container: "docker://sethnr/pgcoe_bacseq:0.01"
    shell:
        """
        samtools view -b -F 4 -F 2048 -o {output.aligned} {input.insorted}
        """

rule sam_subsample:
    input:
        aligned = 'results/align/{sample}_{target}_aligned.bam'
    output:
        subsamp = temporary('results/align/{sample}_{target}_thin.bam'),
        subfactor = 'results/align/{sample}_{target}.substat'
    resources:
        mem_mb=16000,
        runtime=180,
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
        aligned = 'results/align/{sample}_{target}_thin.bam'
    output:
        sorted = 'results/align/{sample}_{target}.bam',
        idx = 'results/align/{sample}_{target}.bam.bai'
    resources:
        mem_mb=16000,
        runtime=180,
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




#rule indexbam:
#    input:
#        bam = '{samplename}.bam'
#    output:
#        indexedbam = '{samplename}.bam.bai'
#    log:
#        stderr="logs/bamindex/{samplename}.err"
#    shell:
#        """
#        samtools index {input.bam} 2>&1 >  {log.stderr}
#        """




### STATS ###

rule depth:
    input:
        bams='results/align/{sample}_{target}.bam',
        subfactor = 'results/align/{sample}_{target}.substat'
    output:
        depth=temporary('results/align/{sample}_{target}_depth.txt'),
        dwins=temporary('results/align/{sample}_{target}_depthwins.txt'),
        dhist=temporary('results/align/{sample}_{target}_depthhist.txt'),
    resources:
        mem_mb=8000,
        runtime=180,
    group: "statsgroup"
    params:
        maxdepth=0,
        minmapqual=60,
        minbasequal=13,
        winsize=10,
        prefix='results/align/{sample}_{target}'
    log:
        stderr="logs/depth/{sample}_{target}.err"
    shell:
        """
        samtools depth -a -H {input.bams} -o {output.depth} 2>&1 >  {log.stderr}

        python scripts/get_depth_distribution.py -d {output.depth} \
            -s {wildcards.sample} -F `cat {input.subfactor} | cut -f 2,2` \
            -w {params.winsize} -o {params.prefix}

        """


rule alignstats:
    input: 
        subfactor = 'results/align/{sample}_{target}.substat',
        flagstats = 'results/align/{sample}_{target}.flagstats',
        dhist='results/align/{sample}_{target}_depthhist.txt',
    output:
        stats='results/align/{sample}_{target}_alignstats.txt',
    group: "statsgroup"
    params:
        mindepth=10
    run:
        sample = wildcards.sample
        target = wildcards.target

        #get subsampling factor
        with open(input.subfactor, "r") as f:
            l = f.read().split()
            subfact = l[1] 
        f.close()

        #get reads aligned
        reads = -1
        aligned = -1
        paired = -1
        with open(input.flagstats, "r") as f:
            for l in f:
                l = l.strip().split('\t')
                passreads = l[0]
                failreads = l[1]
                stat = l[2]
                if stat == 'total (QC-passed reads + QC-failed reads)': 
                    reads = int(passreads)
                elif stat == "mapped": 
                    aligned = int(passreads)
                elif stat == "properly paired": 
                    paired = int(passreads)
        f.close()

        #get coverage / depth
        goodcov = 0
        cov = 0
        gsize = 0
        dtotal = 0
        with open(input.dhist, "r") as f:
            for l in f:
                l = l.split("\t")
                depth = int(l[1])
                count = int(l[2])
                cdepth = int(l[3])
                dtotal += cdepth * count
                if cdepth >= params.mindepth:
                    gsize += count
                    cov += count
                    goodcov += count
                elif cdepth > 0:
                    gsize += count
                    cov += count
                elif cdepth == 0:
                    gsize += count
        f.close()
        meandepth = round(dtotal/gsize)


        #subprocess.run(["samtools","index",output.subsamp])
        f = open(output.stats, "w")
        print("\t".join(map(str,[sample, target, subfact,
	    reads, aligned, paired,
	    meandepth, goodcov, cov, gsize])),file=f)
        f.close()
    