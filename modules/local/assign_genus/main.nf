process ASSIGN_GENUS {
    tag "assign-genus"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path samplesheet
    path gtdbtk_summary
    path genus_map

    output:
    path "sample_assignments.tsv"                 , emit: assignments
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gtdbtk_arg = gtdbtk_summary ? "--gtdbtk-summary ${gtdbtk_summary}" : ""
    def genus_map_arg = genus_map ? "--genus-map ${genus_map}" : ""
    """
    assign_genus.py \\
        --samplesheet ${samplesheet} \\
        --project-dir ${projectDir} \\
        ${gtdbtk_arg} \\
        ${genus_map_arg} \\
        --output sample_assignments.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assign_genus: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    cat <<'TSV' > sample_assignments.tsv
    sample	fasta	genus	genus_source
    MT1012	tests/data/MT1012.fasta	TestGenus	samplesheet
    MT1013	tests/data/MT1013.fasta	TestGenus	samplesheet
    TSV
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assign_genus: 1.0.0
    END_VERSIONS
    """
}
