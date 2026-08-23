process SPECIES_DECISION {
    tag "species-decision"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path assignments
    path query_ani
    path ddh_results
    path pairs

    output:
    path "species_delimitation_summary.tsv", emit: summary
    path "species_delimitation_report.html", emit: report
    path "versions.yml"                    , emit: versions

    script:
    """
    shopt -s nullglob
    if [ -f "${query_ani}" ]; then
        cp "${query_ani}" query_ani_merged.tsv
    elif [ -d "${query_ani}" ]; then
        merge_tsv.py \\
            --inputs ${query_ani}/*.tsv \\
            --output query_ani_merged.tsv
    else
        merge_tsv.py \\
            --inputs ${query_ani} \\
            --output query_ani_merged.tsv
    fi

    if ls *_vs_*.tsv >/dev/null 2>&1; then
        ddh_files=()
        for f in *_vs_*.tsv; do
            if head -n 1 "\$f" | grep -q 'ddh_formula2'; then
                ddh_files+=( "\$f" )
            fi
        done
        if [ \${#ddh_files[@]} -gt 0 ]; then
            merge_tsv.py \\
                --inputs "\${ddh_files[@]}" \\
                --output ddh_merged.tsv
        fi
    elif [ -f "${ddh_results}" ]; then
        cp "${ddh_results}" ddh_merged.tsv
    elif [ -d "${ddh_results}" ]; then
        merge_tsv.py \\
            --inputs ${ddh_results}/*.tsv \\
            --output ddh_merged.tsv
    else
        merge_tsv.py \\
            --inputs ${ddh_results} \\
            --output ddh_merged.tsv
    fi

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        species_decision: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    printf 'sample\\tpredicted_genus\\tclosest_type_strain\\tclosest_type_strain_name\\ttop_fastani_identity\\ttop_ddh_formula2\\tddh_method_used\\ttaxonomic_verdict\\tis_novel_candidate\\nS1\\tGenus\\tGCF.1\\tOrganism\\t96.00\\t72.00\\tlocal\\tnovel_species_candidate\\ttrue\\n' > species_delimitation_summary.tsv
    echo '<html><body>summary</body></html>' > species_delimitation_report.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        species_decision: 1.0.0
    END_VERSIONS
    """
}
