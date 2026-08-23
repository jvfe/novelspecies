process GGDC_INGEST {
    tag "ggdc-ingest"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path results_dir
    path pairs

    output:
    path "ggdc_ddh.tsv", emit: ddh
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ggdc_ingest.py \\
        --results-dir ${results_dir} \\
        --pairs-tsv ${pairs} \\
        --output ggdc_ddh.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ggdc_ingest: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'sample\\tgenus\\tref_accession\\tref_organism\\tddh_formula2\\tdistance_d4\\tprobability_same_species\\tddh_method\\tresult_file\\nS1\\tGenus\\tGCF.1\\tOrganism\\t72.00\\t0.04\\t85.0\\tggdc_web\\tNA\\n' > ggdc_ddh.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ggdc_ingest: 1.0.0
    END_VERSIONS
    """
}
