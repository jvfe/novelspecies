process NCBI_DATASETS_DOWNLOAD {
    tag "$meta.genus"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Slim image: datasets + dataformat + unzip only. No Python.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/staphb/ncbi-datasets:18.35.0' :
        'quay.io/staphb/ncbi-datasets:18.35.0' }"

    input:
    tuple val(meta), val(genus)

    output:
    tuple val(meta), path("ncbi_dataset.zip"), path("assembly_report.tsv"), path("accessions.txt"), emit: dataset
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def api_key_arg   = params.ncbi_api_key ? "--api-key ${params.ncbi_api_key}" : ""
    def reference_arg = params.reference_only ? "--reference" : ""
    def atypical_arg  = params.exclude_atypical ? "--exclude-atypical" : ""
    def annotated_arg = params.annotated_only ? "--annotated" : ""
    def type_only     = params.type_material_only ? '1' : '0'
    def fallback      = params.allow_non_type_fallback ? '1' : '0'
    """
    datasets summary genome taxon "${genus}" \\
        --assembly-source ${params.assembly_source} \\
        --assembly-level ${params.assembly_level} \\
        ${reference_arg} \\
        ${atypical_arg} \\
        ${annotated_arg} \\
        ${api_key_arg} \\
        --as-json-lines \\
        > assembly_summary.jsonl

    if [ ! -s assembly_summary.jsonl ]; then
        echo "ERROR: NCBI returned no assemblies for taxon '${genus}'" >&2
        exit 1
    fi

    dataformat tsv genome \\
        --inputfile assembly_summary.jsonl \\
        --force \\
        --fields accession,organism-name,organism-infraspecific-strain,assminfo-level,assminfo-name,assmstats-total-sequence-len,source_database,type_material-label,type_material-display_text \\
        > assembly_report.tsv

    awk -F '\\t' -v type_only=${type_only} -v fallback=${fallback} '
        function trim(s) { gsub(/\\r/, "", s); return s }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                h = tolower(trim(\$i))
                if (h == "accession" || h == "assembly accession") acc = i
                if (h ~ /type material label/) lab = i
                if (h ~ /type material display/) disp = i
            }
            if (!acc) acc = 1
            next
        }
        {
            a = trim(\$acc)
            if (a == "") next
            print a > "all_accessions.txt"
            is_type = 0
            if (lab) {
                l = toupper(trim(\$lab))
                if (l == "TYPE_MATERIAL") is_type = 1
            }
            if (disp && index(tolower(trim(\$disp)), "type material") > 0) is_type = 1
            if (is_type) print a > "typed_accessions.txt"
        }
    ' assembly_report.tsv

    if [ -s typed_accessions.txt ]; then
        head -n ${params.max_references_per_genus} typed_accessions.txt > accessions.txt
    elif [ "${fallback}" = "1" ] && [ -s all_accessions.txt ]; then
        echo "WARN: no type-material assemblies for ${genus}; using fallback accessions" >&2
        head -n ${params.max_references_per_genus} all_accessions.txt > accessions.txt
    else
        echo "ERROR: no type-material assemblies for taxon '${genus}'" >&2
        exit 1
    fi

    datasets download genome accession \\
        --inputfile accessions.txt \\
        --include genome \\
        ${api_key_arg} \\
        --no-progressbar \\
        --filename ncbi_dataset.zip

    printf '%s\\n' '"${task.process}":' '    ncbi-datasets-cli: "18.35.0"' > versions.yml
    """

    stub:
    """
    echo 'accession' > assembly_report.tsv
    echo 'GCF_STUB.1' >> assembly_report.tsv
    echo 'GCF_STUB.1' > accessions.txt
    echo PK > ncbi_dataset.zip
    printf '%s\\n' '"${task.process}":' '    ncbi-datasets-cli: "18.35.0"' > versions.yml
    """
}
