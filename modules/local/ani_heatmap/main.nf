process ANI_HEATMAP {
    tag "$meta.genus"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyani:0.2.13--pyhdfd78af_0' :
        'quay.io/biocontainers/pyani:0.2.13--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ani_tsv)

    output:
    tuple val(meta), path("*.png") , emit: png
    tuple val(meta), path("*.tsv") , emit: tables
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.genus}"
    """
    export MPLCONFIGDIR="\$PWD/mplconfig"
    mkdir -p "\$MPLCONFIGDIR"

    ani_heatmap.py \\
        --ani-tsv ${ani_tsv} \\
        --genus "${meta.genus}" \\
        --prefix "${prefix}" \\
        --cluster-threshold ${params.ani_cluster_threshold} \\
        --versions versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.genus}"
    """
    touch ${prefix}.ani_heatmap.png ${prefix}.ani_matrix.tsv ${prefix}.ani_linkage.tsv ${prefix}.ani_clusters_${params.ani_cluster_threshold.intValue()}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ani_heatmap: 1.0.0
    END_VERSIONS
    """
}
