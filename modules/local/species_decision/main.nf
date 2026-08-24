process SPECIES_DECISION {
    tag "species-decision"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path assignments
    path query_ani, stageAs: 'ani_inputs/?'
    path ddh_results, stageAs: 'ddh_inputs/?'
    path pairs

    output:
    path "species_delimitation_summary.tsv", emit: summary
    path "species_delimitation_report.html", emit: report
    path "versions.yml"                    , emit: versions

    script:
    """
    shopt -s nullglob
    merge_tsv.py \\
        --inputs ani_inputs/* \\
        --output query_ani_merged.tsv

    if [ ! -s query_ani_merged.tsv ]; then
        printf 'query\\treference\\tani_percent\\tsample\\tgenus\\tref_accession\\n' > query_ani_merged.tsv
    fi

    merge_tsv.py \\
        --inputs ddh_inputs/* \\
        --output ddh_merged.tsv

    if [ ! -s ddh_merged.tsv ]; then
        printf 'sample\\tgenus\\tref_accession\\tref_organism\\tddh_formula2\\tddh_method\\n' > ddh_merged.tsv
    fi

    species_decision.py \\
        --assignments-tsv ${assignments} \\
        --query-ani-tsv query_ani_merged.tsv \\
        --ddh-tsv ddh_merged.tsv \\
        --pairs-tsv ${pairs} \\
        --summary-tsv species_delimitation_summary.tsv \\
        --report-html species_delimitation_report.html \\
        --ani-threshold ${params.ani_threshold} \\
        --ddh-threshold ${params.ddh_threshold} \\
        --borderline-low ${params.ani_borderline_low} \\
        --borderline-high ${params.ani_borderline_high} \\
        --ddh-method ${params.ddh_method}

    printf '%s\\n' '"${task.process}":' '    species_decision: "1.0.0"' > versions.yml
    """

    stub:
    """
    printf 'sample\\tpredicted_genus\\tclosest_type_strain\\tclosest_type_strain_name\\ttop_fastani_identity\\ttop_ddh_formula2\\tddh_method_used\\ttaxonomic_verdict\\tis_novel_candidate\\nS1\\tGenus\\tGCF.1\\tOrganism\\t96.00\\t72.00\\tlocal\\tnovel_species_candidate\\ttrue\\n' > species_delimitation_summary.tsv
    echo '<html><body>summary</body></html>' > species_delimitation_report.html
    printf '%s\\n' '"${task.process}":' '    species_decision: "1.0.0"' > versions.yml
    """
}
