process PREPARE_REFERENCE_LIST {
    tag "$meta.genus"
    label 'process_single'

    input:
    tuple val(meta), val(rows)

    output:
    tuple val(meta), path("references.tsv"), emit: references
    tuple val(meta), path("ref_list.txt") , emit: ref_list
    tuple val(meta), path("staged_refs")  , emit: staged_refs

    script:
    """
    mkdir -p staged_refs
    python3 <<'PY'
import csv
from pathlib import Path
import shutil

rows = ${groovy.json.JsonOutput.toJson(rows)}
fieldnames = [
    "genus", "accession", "organism_name", "strain",
    "assembly_level", "genome_size", "is_type_strain",
    "selection_mode", "fasta",
]
out_rows = []
for row in rows:
    src = Path(row["fasta"])
    dst = Path("staged_refs") / f"{row['accession']}.fna"
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    row = dict(row)
    row["fasta"] = str(dst.resolve())
    out_rows.append(row)

with open("references.tsv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\\t")
    writer.writeheader()
    writer.writerows(out_rows)

with open("ref_list.txt", "w") as handle:
    for row in out_rows:
        handle.write(f"{row['fasta']}\\n")
PY
    """

    stub:
    """
    mkdir staged_refs
    printf 'genus\\taccession\\torganism_name\\tstrain\\tassembly_level\\tgenome_size\\tis_type_strain\\tselection_mode\\tfasta\\n${meta.genus}\\tGCF_STUB.1\\tOrganism\\tTYPE\\tComplete\\t5000000\\ttrue\\treference_sheet\\tstaged_refs/GCF_STUB.1.fna\\n' > references.tsv
    echo 'staged_refs/GCF_STUB.1.fna' > ref_list.txt
    touch staged_refs/GCF_STUB.1.fna
    """
}
