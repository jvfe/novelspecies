process GBDP_DDH {
    tag "$meta.id vs $meta.ref_accession"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyani:0.2.13--pyhdfd78af_0' :
        'quay.io/biocontainers/pyani:0.2.13--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(query), path(reference)

    output:
    tuple val(meta), path("*.tsv"), emit: ddh
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_vs_${meta.ref_accession}"
    def blast_args = task.ext.args ?: '-task blastn -word_size 38 -dust no -evalue 1e-2 -max_target_seqs 1000000'
    """
    gbdp_ddh.py \\
        --query ${query} \\
        --reference ${reference} \\
        --sample ${meta.id} \\
        --genus ${meta.genus} \\
        --ref-accession ${meta.ref_accession} \\
        --ref-organism "${meta.ref_organism ?: ''}" \\
        --blast-args "${blast_args}" \\
        --output ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | head -n 1 | sed 's/.*: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_vs_${meta.ref_accession}"
    """
    printf 'sample\\tgenus\\tref_accession\\tref_organism\\tformula\\thsp_count\\thsp_identity_fraction\\tgbdp_distance_d4\\tddh_formula2\\thsp_coverage_percent\\tddh_method\\n${meta.id}\\t${meta.genus}\\t${meta.ref_accession}\\tNA\\tGGDC_formula_2_d4\\t100\\t0.960000\\t0.040000\\t72.00\\t96.00\\tlocal_blast\\n' > ${prefix}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gbdp_ddh: 1.0.0
    END_VERSIONS
    """
}
