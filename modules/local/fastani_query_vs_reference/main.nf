process FASTANI_QUERY_VS_REFERENCE {
    tag "$meta.genus"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.34--hb66fcc3_7' :
        'quay.io/biocontainers/fastani:1.34--hb66fcc3_7' }"

    input:
    tuple val(meta), path(query_list), path(ref_list), path(references)

    output:
    tuple val(meta), path("*.tsv"), emit: ani
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.genus}_query_vs_reference"
    """
    fastANI \\
        --ql ${query_list} \\
        --rl ${ref_list} \\
        --threads ${task.cpus} \\
        -o ${prefix}.raw.txt

    awk -v comparison="query_vs_reference" 'BEGIN { OFS="\\t"; print "query","reference","ani_percent","comparison_type" }
    NF >= 3 {
        n = split(\$1, qp, "/"); q = qp[n]
        sub(/_genomic\\.fna\$/, "", q); sub(/\\.fna\$/, "", q); sub(/\\.fasta\$/, "", q)
        n = split(\$2, rp, "/"); r = rp[n]
        sub(/_genomic\\.fna\$/, "", r); sub(/\\.fna\$/, "", r); sub(/\\.fasta\$/, "", r)
        ani = \$3; if (ani <= 1) { ani = ani * 100 }
        printf "%s\\t%s\\t%.2f\\t%s\\n", q, r, ani, comparison
    }' ${prefix}.raw.txt > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.genus}_query_vs_reference"
    """
    printf 'query\\treference\\tani_percent\\tcomparison_type\\nS1\\tGCF_REF.1\\t96.50\\tquery_vs_reference\\n' > ${prefix}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: 1.34
    END_VERSIONS
    """
}
