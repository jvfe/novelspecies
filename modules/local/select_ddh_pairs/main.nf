process SELECT_DDH_PAIRS {
    tag "select-ddh-pairs"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path ani_results
    path references
    path assignments

    output:
    path "ddh_pairs.tsv", emit: pairs
    path "versions.yml" , emit: versions

    script:
    """
    shopt -s nullglob
    if [ -f "${ani_results}" ]; then
        cp "${ani_results}" query_ani_merged.tsv
    elif [ -d "${ani_results}" ]; then
        merge_tsv.py \\
            --inputs ${ani_results}/*.tsv \\
            --output query_ani_merged.tsv
    else
        merge_tsv.py \\
            --inputs ${ani_results} \\
            --output query_ani_merged.tsv
    fi

    if [ -f "${references}" ]; then
        cp "${references}" references_merged.tsv
    elif [ -d "${references}" ]; then
        merge_tsv.py \\
            --inputs ${references}/*.tsv \\
            --output references_merged.tsv
    else
        merge_tsv.py \\
            --inputs ${references} \\
            --output references_merged.tsv
    fi

    select_ddh_pairs.py \\
        --ani-tsv query_ani_merged.tsv \\
        --references-tsv references_merged.tsv \\
        --assignments-tsv ${assignments} \\
        --top-n ${params.ddh_top_n} \\
        --output ddh_pairs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        select_ddh_pairs: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'sample\\tgenus\\tquery_fasta\\tref_accession\\tref_organism\\treference_fasta\\tani_percent\\trank\\nS1\\tGenus\\tq.fna\\tGCF.1\\tOrganism\\tr.fna\\t96.0\\t1\\n' > ddh_pairs.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        select_ddh_pairs: 1.0.0
    END_VERSIONS
    """
}
