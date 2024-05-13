configfile: "config/config.yaml"

rule all:
    input:
        nextstrain=expand('results/nextstrain/{target}_auspice.json',target=["RSVA","RSVB"])


rule next_index:
    input:
        fastas='results/summary/{target}_consensus.fasta'
    output:
        index='results/nextstrain/{target}_seqindex.tsv',
    resources:
        mem_mb=4000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        err="logs/nextstrain/index_{target}.err"
    shell:
        """
        augur index \
            --sequences {input.fastas} \
            --output {output.index} >> {log.err}
        """

rule next_filter:
    input:
        fastas='results/summary/{target}_consensus.fasta',
        index='results/nextstrain/{target}_seqindex.tsv',
	    kill='results/nextstrain/{target}_kill-list.txt'
    output:
        filtered='results/nextstrain/{target}_filtered.fasta',
    resources:
        mem_mb=4000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stdout="logs/nextstrain/filter_{target}.out",
        stderr="logs/nextstrain/filter_{target}.err"
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


rule next_align:
    input:
        filtered='results/nextstrain/{target}_filtered.fasta',
        index='results/nextstrain/{target}_seqindex.tsv',
    output:
        aligned='results/nextstrain/{target}_aligned.fasta',
    resources:
        mem_mb=8000,
        runtime=600,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stderr="logs/nextstrain/filter_{target}.err"
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
        aligned='results/nextstrain/{target}_aligned.fasta',
        index='results/nextstrain/{target}_seqindex.tsv',
    output:
        tree='results/nextstrain/{target}_raw.nwk',
    resources:
        mem_mb=8000,
        runtime=240,
    container: "docker://nextstrain/base:build-20240501T170101Z"
    log:
        stderr="logs/nextstrain/tree_{target}.err"
    shell:
        """
        augur tree \
        --alignment {input.aligned} \
        --output {output.tree}
        """


rule next_refine:
    input:
        tree='results/nextstrain/{target}_raw.nwk',
        aligned='results/nextstrain/{target}_aligned.fasta',
    output:
        tree = "results/nextstrain/{target}_tree.nwk",
        node_data = "results/nextstrain/{target}_branchlens.json"
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
        tree = 'results/nextstrain/{target}_tree.nwk',
        alignment = 'results/nextstrain/{target}_aligned.fasta',
    output:
        node_data = "results/nextstrain/{target}_ntmuts.json"
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
        tree =  "results/nextstrain/{target}_tree.nwk",
        branch_lengths = "results/nextstrain/{target}_branchlens.json",
        nt_muts = "results/nextstrain/{target}_ntmuts.json",
    output:
        auspice_json = "results/nextstrain/{target}_auspice.json",
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

