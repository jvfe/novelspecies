process NCBI_DATASETS_STAGE {
    tag "$meta.genus"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(dataset_zip), path(assembly_report)

    output:
    tuple val(meta), path("references.tsv"), emit: references
    tuple val(meta), path("ref_list.txt") , emit: ref_list
    tuple val(meta), path("staged_refs")  , emit: staged_refs
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def type_only = params.type_material_only ? '--type-material-only' : ''
    def fallback  = params.allow_non_type_fallback ? '--allow-non-type-fallback' : ''
    """
    python3 - <<'PY'
import zipfile
from pathlib import Path
archive = Path("${dataset_zip}")
if not archive.exists() or archive.stat().st_size == 0:
    raise SystemExit("ERROR: NCBI dataset zip is missing or empty")
with zipfile.ZipFile(archive) as handle:
    handle.extractall()
PY

    prepare_reference_manifest.py \\
        --report-tsv ${assembly_report} \\
        --data-dir ncbi_dataset/data \\
        --genus "${meta.genus}" \\
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
        ncbi_datasets_stage: 1.0.0
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
        ncbi_datasets_stage: 1.0.0
    END_VERSIONS
    """
}
