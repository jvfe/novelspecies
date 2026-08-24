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
    def from_type_arg = params.type_material_only ? "--from-type" : ""
    def type_only     = params.type_material_only ? '1' : '0'
    def fallback      = params.allow_non_type_fallback ? '1' : '0'
    def max_refs      = params.max_references_per_genus
    """
    genus="${genus}"
    # GTDB placeholder genera (WGS project IDs) are not NCBI taxa, e.g. CAJYPV01, JAADFP01
    if echo "\$genus" | grep -Eq '^[A-Z]{2,10}[0-9]{2,}\$|^UBA[0-9]+\$|^GCA-[0-9]+\$'; then
        echo "SKIP: '\$genus' is a GTDB placeholder genus (WGS/MAG identifier), not an NCBI taxon." >&2
        echo "      Provide a --reference_sheet for this lineage, or drop those queries." >&2
        exit 4
    fi

    # GTDB satellite genera: Enterococcus_B, Bacillus_AE, Paenibacillus_O
    stripped=\$(echo "\$genus" | sed -E 's/_[A-Z]+\$//')

    : > summary.err

    pick_prokaryote_taxid() {
        # datasets prints: Bacillus (genus, taxid: 55087, walking sticks)
        #                  Bacillus (genus, taxid: 1386, firmicutes)
        grep -oE 'taxid: [0-9]+, [^)]+' "\$1" \\
            | grep -iE 'firmicutes|proteobacteria|bacter|actinobacteria|spirochetes|cyanobacteria|chlamydiae|tenericutes|fusobacteria|deinococci|thermotogae|planctomycetes|verrucomicrobia|acidobacteria|chloroflexi|synergistetes' \\
            | grep -oE '[0-9]+' \\
            | head -n 1
    }

    run_datasets_summary() {
        taxon="\$1"
        shift
        rm -f assembly_summary.jsonl summary.try.err
        datasets summary genome taxon "\$taxon" \\
            --assembly-source ${params.assembly_source} \\
            --limit ${max_refs} \\
            --mag exclude \\
            ${from_type_arg} \\
            ${reference_arg} \\
            ${atypical_arg} \\
            ${annotated_arg} \\
            ${api_key_arg} \\
            --as-json-lines \\
            "\$@" \\
            > assembly_summary.jsonl 2> summary.try.err
        cat summary.try.err >> summary.err
        [ -s assembly_summary.jsonl ]
    }

    try_summary() {
        taxon="\$1"
        shift
        echo "NCBI summary: taxon='\$taxon' \$*" >&2
        if run_datasets_summary "\$taxon" "\$@"; then
            return 0
        fi
        if grep -q "more than one taxid" summary.try.err; then
            taxid=\$(pick_prokaryote_taxid summary.try.err)
            if [ -n "\$taxid" ]; then
                echo "WARN: '\$taxon' matches multiple NCBI taxids; using bacterial taxid \$taxid" >&2
                run_datasets_summary "\$taxid" "\$@"
                return \$?
            fi
        fi
        return 1
    }

    resolve_genus=""
    for taxon in "\$genus" "\$stripped"; do
        if [ "\$taxon" = "\$resolve_genus" ] || [ -z "\$taxon" ]; then
            continue
        fi
        if try_summary "\$taxon" --assembly-level ${params.assembly_level}; then
            resolve_genus="\$taxon"
            break
        fi
        # Draft-only type strains (e.g. Mangrovactinospora gilvigrisea) are contig-level
        if try_summary "\$taxon"; then
            echo "WARN: no ${params.assembly_level} assemblies for '\$taxon'; including contig-level genomes" >&2
            resolve_genus="\$taxon"
            break
        fi
        resolve_genus="\$taxon"
    done

    if [ ! -s assembly_summary.jsonl ] && [ "${params.reference_only}" = "false" ]; then
        echo "WARN: type-material summary failed; retrying with NCBI reference genomes only" >&2
        for taxon in "\$genus" "\$stripped"; do
            if [ -z "\$taxon" ]; then
                continue
            fi
            if try_summary "\$taxon" --reference --assembly-level ${params.assembly_level}; then
                echo "WARN: NCBI has no type-material package for '\$genus'; using reference genomes of '\$taxon'" >&2
                resolve_genus="\$taxon"
                break
            fi
        done
    fi

    if [ ! -s assembly_summary.jsonl ]; then
        echo "SKIP: NCBI has no ${params.assembly_source} assemblies for taxon '\$genus'." >&2
        cat summary.err >&2
        exit 4
    fi

    if [ "\$resolve_genus" != "\$genus" ]; then
        echo "WARN: NCBI has no assemblies for GTDB genus '\$genus'; using '\$resolve_genus'" >&2
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
                if (index(l, "TYPE_MATERIAL") > 0) is_type = 1
            }
            if (disp && index(tolower(trim(\$disp)), "type material") > 0) is_type = 1
            if (is_type) print a > "typed_accessions.txt"
        }
    ' assembly_report.tsv

    # --from-type already restricted the summary; keep those accessions even if
    # dataformat left the type-material columns empty.
    if [ ! -s typed_accessions.txt ] && [ -s all_accessions.txt ] && [ "${type_only}" = "1" ]; then
        cp all_accessions.txt typed_accessions.txt
    fi

    if [ -s typed_accessions.txt ]; then
        head -n ${max_refs} typed_accessions.txt > accessions.txt
    elif [ "${fallback}" = "1" ] && [ -s all_accessions.txt ]; then
        echo "WARN: no type-material assemblies for ${genus}; using fallback accessions" >&2
        head -n ${max_refs} all_accessions.txt > accessions.txt
    else
        echo "SKIP: no type-material assemblies for taxon '${genus}'" >&2
        exit 4
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
