# rule all:
#     input:
#         auspice_json = "auspice/zika.json",
#         input_fasta = "data/sequences.fasta",
#         input_metadata = "data/metadata.tsv",
#         dropped_strains = "config/dropped_strains.txt",
#         reference = "config/zika_outgroup.gb",
#         colors = "config/colors.tsv",
#         lat_longs = "config/lat_longs.tsv",
#         auspice_config = "config/auspice_config.json"



rule next_index:
    input:
        fastas='results/summary/{target}-consensus.fasta'
    output:
        index='results/nextstrain/{target}-seqindex.tsv',
    resources:
        mem_mb=4000,
        runtime=240,
    container: docker://nextstrain/base:build-20230525T143814Z
    log:
        stderr="logs/nextstrain/index-{target}.err"
    shell:
        """
        augur index \
            --sequences {input.fastas} \
            --output {output.index} >> {log.err}
        """

rule next_filter:
    input:
        fastas='results/ivar/{target}-consensus.fasta',
        index='results/nextstrain/{target}-seqindex.tsv',
    output:
        filtered='results/nextstrain/{target}-filtered.fasta',
    resources:
        mem_mb=4000,
        runtime=240,
    container: docker://nextstrain/base:build-20230525T143814Z
    log:
        stderr="logs/nextstrain/filter-{target}.err"
    params:
        mindate=2000,
        seq_per_group=50,
        ref=os.path.join(config['refsdir'],wildcards.target,".fa"),
        metadata=config['metadata'],
        dropped=config['nextdropped'],
        group_by = "country year month",
    shell:
        """
        augur filter \
            --sequences {input.fastas} \
            --sequence-index {input.index} \
            --metadata {params.metadata} \
            --exclude {params.dropped} \
            --group-by {params.group_by} \
            --sequences-per-group {params.seq_per_group} \
            --min-date {params.mindate} \
            --output {output.filtered}
        """


rule next_align:
    input:
        filtered='results/nextstrain/{target}-filtered.fasta',
        index='results/nextstrain/{target}-seqindex.tsv',
    output:
        aligned='results/nextstrain/{target}-aligned.fasta',
    resources:
        mem_mb=4000,
        runtime=240,
    container: docker://nextstrain/base:build-20230525T143814Z
    log:
        stderr="logs/nextstrain/filter-{target}.err"
    params:
        mindate=2000,
        seq_per_group=50,
        ref=os.path.join(config['refsdir'],wildcards.target,".fa")
    shell:
        """
        augur align \
            --sequences {output.filtered} \
            --reference-sequence {params.ref} \
            --output {output.aligned} \
            --fill-gaps
        """


rule next_tree:
    input:
        aligned='results/nextstrain/{target}-aligned.fasta',
        index='results/nextstrain/{target}-seqindex.tsv',
    output:
        raw='results/nextstrain/{target}-raw.nwk',
    resources:
        mem_mb=8000,
        runtime=240,
    container: docker://nextstrain/base:build-20230525T143814Z
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
        raw='results/nextstrain/{target}-raw.nwk',
        aligned='results/nextstrain/{target}-aligned.fasta',
        metadata = input_metadata
    output:
        tree = "results/nextstrain/{target}-tree.nwk",
        node_data = "results/nextstrain/{target}-branchlens.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    resources:
        mem_mb=8000,
        runtime=240,
    container: docker://nextstrain/base:build-20230525T143814Z
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule next_ancestral:
    input:
        tree = 'results/nextstrain/{target}-tree.nwk',
        alignment = 'results/nextstrain/{target}-alignment.fasta',
    output:
        node_data = "results/nextstrain/{target}-ntmuts.json"
    params:
        inference = "joint"
    resources:
        mem_mb=8000,
        runtime=240,
    container: docker://nextstrain/base:build-20230525T143814Z
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
    container: docker://nextstrain/base:build-20230525T143814Z
    params:
        metadata=config['metadata'],
        colors = config['nscolors'],
        lat_longs = config['nslatlongs'],
        auspice_config = config['auspice_config']
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {params.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --colors {params.colors} \
            --lat-longs {params.lat_longs} \
            --auspice-config {params.auspice_config} \
            --include-root-sequence \
            --output {output.auspice_json}
        """

