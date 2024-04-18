rule run_mash:
    input:
        msh=config['refsdir']+'/all.msh'


rule mash_index:
    input:
        fastas=expand(config['refsdir']+"/{target}.fasta",target=TARGETS),
    output:
        msh=config['refsdir']+"/all.msh"
    log:
        log="logs/mash/all.mash.log"
    params:
        idx_script = "scripts/indexer.sh",
        genome_size = "11k",
        refdir=config["refsdir"],
        localfastas=expand("{target}.fasta",target=TARGETS)
    resources:
        partition="day",
        mem_mb="40G",
        cpus_per_task=4,
        runtime=300
    container: "docker://sethnr/pgcoe_anypipe:0.01"
    shell:"""
        cd {input.refdir};
        {params.idx_script} \
                -m {output.msh} -g {params.genome_size} \
                {params.localfastas} >> {log.log} 2>&1
        """

rule mash_calltarget:
    input:
        msh = expand("{refsdir}/all.msh",refsdir=config["refsdir"]),
        read_location = "{readsdir}/{sample}"
    output:
        mashout = "results/mash/{sample}_mash.txt",
        mashcalls = "results/mash/{sample}_calls.txt"
    resources:
        partition="day",
        mem_mb="8G",
        cpus_per_task=1,
        runtime=300
    container: "docker://sethnr/pgcoe_anypipe:0.01"
    params:
        reads=10000, # compare top N reads to refs
        bloom=10,    # bloom filter kmers with < N coverage (seq errors)
        gsize="11k", # estimated genome size (for prob assignment)
        prob=1e-50,  # max mash prob to call
        dist=0.25,   # max mash dist to call
        masher = "scripts/masher.sh",
        prefix="results/mash/{sample}",
        refdir=config["refsdir"],
    log:
        "logs/getstrain_{sample}.log",
    shell:
        """
        {params.masher} -f {params.refdir}/all.msh -s {wildcards.sample}\
                -r {params.reads} -b {params.bloom} -g {params.gsize} \
                -d {params.dist} -p {params.prob} \
        -o {params.prefix} \
        {input.read_location}
        """

checkpoint: mash_merge_calls:
    input:
        expand("results/mash/{sample}_calls.txt",sample=SAMPLES)
    output:
        "results/mash/all_calls.txt",
    log:
        "logs/mash/all_calls.err"
    shell:
    """ cat {input} 1> {output} 2> {log}"""



# def get_mash_targets(wildcards):
#     with checkpoints.mash_calltarget.get(sample=wildcards.sample,outdir=wildcards.outdir).output.mashcalls.open() as f:
#         mytargets = f.read().splitlines()
#     return mytargets

def get_mash_targets(wildcards):
    with checkpoints.mash_merge_calls.get().output.mashcalls.open() as f:
        mytargets = [L.split[1] for L in f.read().splitlines() if L.split[1] eq wildcards.sample]
    return mytargets

def get_mash_samples(wildcards):
    with checkpoints.mash_merge_calls.get().output.mashcalls.open() as f:
        mysamples = [L.split[0] for L in f.read().splitlines() if L.split[1] eq wildcards.target]
    return mysamples

def get_valid_targets():
    with checkpoints.mash_merge_calls.get().output.mashcalls.open() as f:
        alltargets = list(set([L.split[1] for L in f.read().splitlines()]))
    return alltargets


def get_mash_bams(wildcards):
    mytargets = get_mash_targets(wildcards)
    return expand("results/bams/{sample}_{target}.bam",sample=wildcards.sample,target=mytargets)


