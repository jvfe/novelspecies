process PREPARE_QUERY_LIST {
    tag "$meta.genus"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(assignments)

    output:
    tuple val(meta), path("query_list.txt"), emit: query_list

    script:
    """
    prepare_query_list.py \\
        --assignments ${assignments} \\
        --genus "${meta.genus}" \\
        --output query_list.txt
    """

    stub:
    """
    mkdir -p queries
    echo "queries/query.fna" > query_list.txt
    touch queries/query.fna
    """
}
