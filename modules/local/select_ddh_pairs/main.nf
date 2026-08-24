process SELECT_DDH_PAIRS {
    tag "select-ddh-pairs"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path ani_results, stageAs: 'ani_inputs/?'
    path references, stageAs: 'ref_inputs/?'
    path assignments

    output:
    path "ddh_pairs.tsv", emit: pairs
    path "versions.yml" , emit: versions

    script:
    """
    shopt -s nullglob
    merge_tsv.py \\
        --inputs ani_inputs/* \\
        --output query_ani_merged.tsv

    merge_tsv.py \\
        --inputs ref_inputs/* \\
        --output references_merged.tsv

    select_ddh_pairs.py \\
        --ani-tsv query_ani_merged.tsv \\
        --references-tsv references_merged.tsv \\
        --assignments-tsv ${assignments} \\
        --top-n ${params.ddh_top_n} \\
        --output ddh_pairs.tsv

    printf '%s\\n' '"${task.process}":' '    select_ddh_pairs: "1.0.0"' > versions.yml
    """

    stub:
    """
    printf 'sample\\tgenus\\tquery_fasta\\tref_accession\\tref_organism\\treference_fasta\\tani_percent\\trank\\nS1\\tGenus\\tq.fna\\tGCF.1\\tOrganism\\tr.fna\\t96.0\\t1\\n' > ddh_pairs.tsv
    printf '%s\\n' '"${task.process}":' '    select_ddh_pairs: "1.0.0"' > versions.yml
    """
}
