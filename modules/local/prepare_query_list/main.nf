process PREPARE_QUERY_LIST {
    tag "$meta.genus"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), val(sample_ids), path(fastas)

    output:
    tuple val(meta), path("query_list.txt"), path("queries"), emit: query_list

    script:
    def samples = sample_ids instanceof List ? sample_ids : [sample_ids]
    def fasta_files = fastas instanceof List ? fastas : [fastas]
    def mapping = [samples, fasta_files].transpose().collect { id, fasta -> "${id}\t${fasta.name}" }.join('\n')
    """
    cat > sample_map.tsv <<'EOF'
sample	fasta_name
${mapping}
EOF

    prepare_query_list.py \\
        --mapping sample_map.tsv \\
        --output query_list.txt
    """

    stub:
    def samples = sample_ids instanceof List ? sample_ids : [sample_ids]
    def first = samples[0]
    """
    mkdir -p queries
    echo "queries/${first}.fna" > query_list.txt
    touch queries/${first}.fna
    """
}
