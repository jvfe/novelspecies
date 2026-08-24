process REFERENCE_SHEET_TO_MANIFEST {
    tag "reference-sheet"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path reference_sheet

    output:
    path "genus_manifests", emit: manifests

    script:
    """
    reference_sheet_to_manifest.py \\
        --reference-sheet ${reference_sheet} \\
        --project-dir ${projectDir} \\
        --outdir genus_manifests
    """

    stub:
    """
    mkdir -p genus_manifests/TestGenus/staged_refs
    printf 'genus\\taccession\\torganism_name\\tstrain\\tassembly_level\\tgenome_size\\tis_type_strain\\tselection_mode\\tfasta\\nTestGenus\\tGCF_TYPE.1\\tTestGenus type strain\\tTYPE\\tComplete Genome\\t5000000\\ttrue\\treference_sheet\\tstaged_refs/GCF_TYPE.1.fna\\n' > genus_manifests/TestGenus/references.tsv
    echo 'staged_refs/GCF_TYPE.1.fna' > genus_manifests/TestGenus/ref_list.txt
    touch genus_manifests/TestGenus/staged_refs/GCF_TYPE.1.fna
    """
}
