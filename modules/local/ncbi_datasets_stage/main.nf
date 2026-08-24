process NCBI_DATASETS_STAGE {
    tag "$meta.genus"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(dataset_zip), path(assembly_report), path(accessions)

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
header = archive.read_bytes()[:4]
if header[:2] != b"PK":
    preview = archive.read_text(errors="replace")[:500]
    raise SystemExit(f"ERROR: NCBI download is not a zip file: {preview}")
with zipfile.ZipFile(archive) as handle:
    handle.extractall()
PY

    data_dir=ncbi_dataset/data
    if [ ! -d "\$data_dir" ]; then
        data_dir=.
    fi

    prepare_reference_manifest.py \\
        --report-tsv ${assembly_report} \\
        --data-dir "\$data_dir" \\
        --genus "${meta.genus}" \\
        --downloaded-accessions ${accessions} \\
        ${type_only} \\
        ${fallback} \\
        --max-references ${params.max_references_per_genus} \\
        --staged-dir staged_refs \\
        --references-tsv references.tsv \\
        --ref-list ref_list.txt

    printf '%s\\n' '"${task.process}":' '    ncbi_datasets_stage: "1.0.0"' > versions.yml
    """

    stub:
    """
    mkdir -p staged_refs
    printf 'genus\\taccession\\torganism_name\\tstrain\\tassembly_level\\tgenome_size\\tis_type_strain\\tselection_mode\\tfasta\\n${meta.genus}\\tGCF_STUB.1\\t${meta.genus} sp.\\tTYPE\\tComplete Genome\\t5000000\\ttrue\\ttype_material\\tstaged_refs/GCF_STUB.1.fna\\n' > references.tsv
    echo 'staged_refs/GCF_STUB.1.fna' > ref_list.txt
    touch staged_refs/GCF_STUB.1.fna
    printf '%s\\n' '"${task.process}":' '    ncbi_datasets_stage: "1.0.0"' > versions.yml
    """
}
