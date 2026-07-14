configfile: "config/config.yaml"

rule nextstrain_all:
    input:
        nextstrain=expand('results/nextstrain/{target}-auspice.json',target=["RSVA","RSVB"])


rule next_index:
    input:
        fastas='{fastafile}.fasta'
    output:
        index='{fastafile}-index.tsv',
    resources:
        mem_mb=4000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        err="logs/nextstrain/index-{target}.err"
    shell:
        """
        augur index \
            --sequences {input.fastas} \
            --output {output.index} >> {log.err}
        """

rule next_filter_bgr:
    input:
        fastas='results/nextstrain/{target}-bgr.fasta',
        index='results/nextstrain/{target}-bgr-index.tsv',
	    kill='results/nextstrain/{target}-bgr-kill-list.txt'
    output:
        filtered='results/nextstrain/{target}-bgr-filtered.fasta',
    resources:
        mem_mb=4000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stdout="logs/nextstrain/filter-{target}.out",
        stderr="logs/nextstrain/filter-{target}.err"
    params:
        mindate=2020,
        seq_per_group=10,
        ref=os.path.join(config['refsdir'],"{target}.fa"),
        metadata=config['nsmetadata'],
        group_by = "country year month",
        idcols = "Seq-ID Original-ID name"
    shell:
        """
        augur filter \
            --sequences {input.fastas} \
            --sequence-index {input.index} \
            --metadata {params.metadata} \
            --metadata-id-columns {params.idcols} \
            --exclude {input.kill} \
            --group-by {params.group_by} \
            --sequences-per-group {params.seq_per_group} \
            --min-date {params.mindate} \
            --output {output.filtered} 1> {log.stdout} 2> {log.stderr}
         """




rule next_filter:
    input:
        fastas='results/summary/{target}-consensus.fasta',
        index='results/nextstrain/{target}-consensus-index.tsv',
	    includes='results/nextstrain/{target}-calls.txt'
    output:
        filtered='results/nextstrain/{target}-consensus-filtered.fasta',
    resources:
        mem_mb=4000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stdout="logs/nextstrain/filter-{target}.out",
        stderr="logs/nextstrain/filter-{target}.err"
    params:
        mindate=2020,
        ref=os.path.join(config['refsdir'],"{target}.fa"),
        metadata=config['nsmetadata'],
        idcols = "Seq-ID Original-ID name"
    shell:
        """
        augur filter \
            --sequences {input.fastas} \
            --sequence-index {input.index} \
            --metadata {params.metadata} \
            --include {input.includes} \
            --exclude-all \
            --output {output.filtered} 1> {log.stdout} 2> {log.stderr}
         """

rule cat_fasta:
    input:
        bgr='results/nextstrain/{target}-bgr-filtered.fasta',
        consensus='results/nextstrain/{target}-consensus-filtered.fasta',
    output:
        all='results/nextstrain/{target}-all-filtered.fasta',
    resources:
        mem_mb=8000,
        runtime=600,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stderr="logs/nextstrain/filter-{target}.err"
    params:
        mindate=2000,
        seq_per_group=50,
        ref=os.path.join(config['refsdir'],"{target}.fasta")
    shell:
        """
        cat {input.bgr} {input.consensus} > {output.all}
        """


rule next_align:
    input:
        filtered='results/nextstrain/{target}-bgr-filtered.fasta',
        filtered='results/nextstrain/{target}-filtered.fasta',
        index='results/nextstrain/{target}-seqindex.tsv',
    output:
        aligned='results/nextstrain/{target}-aligned.fasta',
    resources:
        mem_mb=8000,
        runtime=600,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stderr="logs/nextstrain/filter-{target}.err"
    params:
        mindate=2000,
        seq_per_group=50,
        ref=os.path.join(config['refsdir'],"{target}.fasta")
    shell:
        """
        augur align \
            --sequences {input.filtered} \
            --reference-sequence {params.ref} \
            --output {output.aligned} \
            --fill-gaps
        """


rule next_tree:
    input:
        aligned='results/nextstrain/{target}-aligned.fasta',
        index='results/nextstrain/{target}-seqindex.tsv',
    output:
        tree='results/nextstrain/{target}-raw.nwk',
    resources:
        mem_mb=8000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stderr="logs/nextstrain/tree-{target}.err"
    shell:
        """
        augur tree \
        --alignment {input.aligned} \
        --output {output.tree}
        """


rule next_refine:
    input:
        tree='results/nextstrain/{target}-raw.nwk',
        aligned='results/nextstrain/{target}-aligned.fasta',
    output:
        tree = "results/nextstrain/{target}-tree.nwk",
        node_data = "results/nextstrain/{target}-branchlens.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4,
	metadata = config['nsmetadata'],
	idcols = "Seq-ID Original-ID name"
    resources:
        mem_mb=8000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.aligned} \
            --metadata {params.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --metadata-id-columns {params.idcols} \

        """

rule next_ancestral:
    input:
        tree = 'results/nextstrain/{target}-tree.nwk',
        alignment = 'results/nextstrain/{target}-aligned.fasta',
    output:
        node_data = "results/nextstrain/{target}-ntmuts.json"
    params:
        inference = "joint"
    resources:
        mem_mb=8000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """


rule next_export:
    input:
        tree =  "results/nextstrain/{target}-tree.nwk",
        branch_lengths = "results/nextstrain/{target}-branchlens.json",
        nt_muts = "results/nextstrain/{target}-ntmuts.json",
    output:
        auspice_json = "results/nextstrain/{target}-auspice.json",
    resources:
        mem_mb=8000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    params:
        metadata=config['nsmetadata'],
        colors = config['nscolors'],
        lat_longs = config['nslatlongs'],
        traits = config['nstraits'],
        auspice_config = config['nsconfig'],
	idcols = "Seq-ID Original-ID name"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {params.metadata} \
            --node-data {input.branch_lengths} {params.traits} {input.nt_muts} \
            --colors {params.colors} \
            --lat-longs {params.lat_longs} \
            --auspice-config {params.auspice_config} \
            --include-root-sequence \
            --output {output.auspice_json}
        """

