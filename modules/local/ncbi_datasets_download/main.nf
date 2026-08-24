process NCBI_DATASETS_DOWNLOAD {
    tag "$meta.genus"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Galaxy Depot has no working ncbi-datasets-cli 18.x SIF; StaphB image works with Apptainer/Singularity via ociAutoPull
    container 'quay.io/staphb/ncbi-datasets:18.35.0'

    input:
    tuple val(meta), val(genus)

    output:
    tuple val(meta), path("references.tsv"), emit: references
    tuple val(meta), path("ref_list.txt") , emit: ref_list
    tuple val(meta), path("staged_refs")  , emit: staged_refs
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def api_key_arg   = params.ncbi_api_key ? "--api-key ${params.ncbi_api_key}" : ""
    def reference_arg = params.reference_only ? "--reference" : ""
    def atypical_arg  = params.exclude_atypical ? "--exclude-atypical" : ""
    def annotated_arg = params.annotated_only ? "--annotated" : ""
    def type_only     = params.type_material_only ? '--type-material-only' : ''
    def fallback      = params.allow_non_type_fallback ? '--allow-non-type-fallback' : ''
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

    prepare_reference_manifest.py \\
        --report-tsv assembly_report.tsv \\
        --genus "${genus}" \\
        ${type_only} \\
        ${fallback} \\
        --max-references ${params.max_references_per_genus} \\
        --accessions-out accessions.txt

    datasets download genome accession \\
        --inputfile accessions.txt \\
        --include genome \\
        ${api_key_arg} \\
        --no-progressbar \\
        --filename ncbi_dataset.zip

    if command -v python3 >/dev/null 2>&1; then
        python3 - <<'PY'
import zipfile
from pathlib import Path
archive = Path("ncbi_dataset.zip")
if not archive.exists() or archive.stat().st_size == 0:
    raise SystemExit("ERROR: NCBI dataset zip is missing or empty")
with zipfile.ZipFile(archive) as handle:
    handle.extractall()
PY
    elif command -v unzip >/dev/null 2>&1; then
        unzip -qq ncbi_dataset.zip
    else
        echo "ERROR: need python3 or unzip to extract ncbi_dataset.zip" >&2
        exit 127
    fi

    prepare_reference_manifest.py \\
        --report-tsv assembly_report.tsv \\
        --data-dir ncbi_dataset/data \\
        --genus "${genus}" \\
        ${type_only} \\
        ${fallback} \\
        --max-references ${params.max_references_per_genus} \\
        --references-tsv references.raw.tsv \\
        --ref-list ref_list.raw.txt

    mkdir -p staged_refs
    tail -n +2 references.raw.tsv | while IFS=\$'\\t' read -r genus acc organism strain level size is_type mode fasta; do
        cp "\$fasta" "staged_refs/\${acc}.fna"
    done

    python3 <<'PY'
import csv
from pathlib import Path

rows = list(csv.DictReader(open("references.raw.tsv"), delimiter="\\t"))
with open("references.tsv", "w", newline="") as out_handle:
    writer = csv.DictWriter(out_handle, fieldnames=rows[0].keys(), delimiter="\\t")
    writer.writeheader()
    for row in rows:
        staged = Path("staged_refs") / f"{row['accession']}.fna"
        row["fasta"] = staged.as_posix()
        writer.writerow(row)

with open("ref_list.txt", "w") as handle:
    for row in rows:
        handle.write(f"{row['fasta']}\\n")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets-cli: \$(datasets --version 2>&1 | sed 's/datasets version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p staged_refs
    printf 'genus\\taccession\\torganism_name\\tstrain\\tassembly_level\\tgenome_size\\tis_type_strain\\tselection_mode\\tfasta\\n${meta.genus}\\tGCF_STUB.1\\t${meta.genus} sp.\\tTYPE\\tComplete Genome\\t5000000\\ttrue\\ttype_material\\tstaged_refs/GCF_STUB.1.fna\\n' > references.tsv
    echo 'staged_refs/GCF_STUB.1.fna' > ref_list.txt
    touch staged_refs/GCF_STUB.1.fna
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets-cli: 18.35.0
    END_VERSIONS
    """
}
