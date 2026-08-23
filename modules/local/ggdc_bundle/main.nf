process GGDC_BUNDLE {
    tag "ggdc-bundle"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path pairs

    output:
    path "ggdc_manifest.tsv"      , emit: manifest
    path "ggdc_submissions.tar.gz", emit: archive
    path "GGDC_INSTRUCTIONS.txt"  , emit: instructions
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ggdc_bundle.py \\
        --pairs-tsv ${pairs} \\
        --email "${params.ggdc_email ?: ''}" \\
        --bundle-dir ggdc_jobs \\
        --manifest ggdc_manifest.tsv \\
        --archive ggdc_submissions.tar.gz \\
        --instructions GGDC_INSTRUCTIONS.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ggdc_bundle: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch ggdc_manifest.tsv ggdc_submissions.tar.gz GGDC_INSTRUCTIONS.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ggdc_bundle: 1.0.0
    END_VERSIONS
    """
}
