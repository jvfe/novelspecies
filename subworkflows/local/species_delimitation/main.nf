/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPECIES DELIMITATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP                      } from '../../../modules/nf-core/gunzip/main'
include { ASSIGN_GENUS                } from '../../../modules/local/assign_genus/main'
include { NCBI_DATASETS_DOWNLOAD      } from '../../../modules/local/ncbi_datasets_download/main'
include { REFERENCE_SHEET_TO_MANIFEST } from '../../../modules/local/reference_sheet_to_manifest/main'
include { PREPARE_QUERY_LIST          } from '../../../modules/local/prepare_query_list/main'
include { FASTANI_QUERY_VS_REFERENCE  } from '../../../modules/local/fastani_query_vs_reference/main'
include { FASTANI_INTRA_GENUS         } from '../../../modules/local/fastani_intra_genus/main'
include { ANI_HEATMAP                 } from '../../../modules/local/ani_heatmap/main'
include { SELECT_DDH_PAIRS            } from '../../../modules/local/select_ddh_pairs/main'
include { GBDP_DDH                      } from '../../../modules/local/gbdp_ddh/main'
include { GGDC_BUNDLE                   } from '../../../modules/local/ggdc_bundle/main'
include { GGDC_INGEST                   } from '../../../modules/local/ggdc_ingest/main'
include { SPECIES_DECISION              } from '../../../modules/local/species_decision/main'

workflow SPECIES_DELIMITATION {

    take:
    ch_samplesheet    // channel: [ val(meta), path(fasta) ]
    samplesheet_path  // path: original samplesheet CSV
    outdir

    main:
    ch_versions = channel.empty()

    if (params.skip_reference_download && !params.reference_sheet) {
        error("Reference retrieval disabled (--skip_reference_download) but no --reference_sheet was provided.")
    }

    if (params.ddh_method == 'pyani') {
        log.warn "--ddh_method pyani uses the local BLAST GBDP formula-2 engine (not the pyani package)."
    }

    //
    // Decompress gzipped assemblies when required
    //
    ch_input = ch_samplesheet.branch { _meta, fasta ->
        gz:   fasta.name.endsWith('.gz')
        flat: true
    }

    GUNZIP(
        ch_input.gz
            .map { meta, fasta -> [ [ id: meta.id ], fasta ] }
    )

    ch_genomes = ch_input.flat
        .mix(
            GUNZIP.out.gunzip
                .map { meta, fasta ->
                    def sample_meta = [ id: meta.id, single_end: true ]
                    [ sample_meta, fasta ]
                }
        )

    //
    // Assign genus to each query genome
    //
    ch_gtdbtk = params.gtdbtk_summary ? file(params.gtdbtk_summary, checkIfExists: true) : []
    ch_genus_map = params.genus_map ? file(params.genus_map, checkIfExists: true) : []

    ASSIGN_GENUS(
        samplesheet_path,
        ch_gtdbtk,
        ch_genus_map,
    )
    ch_versions = ch_versions.mix(ASSIGN_GENUS.out.versions)

    ch_assigned = ASSIGN_GENUS.out.assignments
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def meta = [
                id: row.sample,
                genus: row.genus,
                genus_source: row.genus_source,
                single_end: true,
            ]
            [ row.sample, meta ]
        }

    ch_query_genomes = ch_genomes
        .map { meta, fasta -> [ meta.id, meta, fasta ] }
        .join(ch_assigned, by: 0)
        .map { _id, meta, fasta, assigned_meta ->
            [ assigned_meta, fasta ]
        }

    ch_by_genus = ch_query_genomes
        .map { meta, fasta -> [ meta.genus, meta, fasta ] }
        .groupTuple(by: 0)
        .map { genus, metas, fastas ->
            def meta = [ id: genus.replaceAll('\\s+', '_'), genus: genus ]
            [ meta, metas, fastas ]
        }

    PREPARE_QUERY_LIST(
        ch_by_genus.map { meta, metas, fastas ->
            [ meta, metas.collect { it.id }, fastas ]
        }
    )

    //
    // Retrieve or load reference genomes per genus
    //
    ch_genus_for_refs = ch_by_genus.map { meta, _metas, _fastas -> [ meta, meta.genus ] }

    ch_reference_sheet_input = params.reference_sheet
        ? channel.value(file(params.reference_sheet, checkIfExists: true))
        : channel.empty()

    REFERENCE_SHEET_TO_MANIFEST(
        ch_reference_sheet_input,
    )

    ch_ncbi_input = params.reference_sheet || params.skip_reference_download
        ? channel.empty()
        : ch_genus_for_refs

    NCBI_DATASETS_DOWNLOAD(
        ch_ncbi_input,
    )
    ch_versions = ch_versions.mix(NCBI_DATASETS_DOWNLOAD.out.versions)

    ch_reference_tables_from_sheet = REFERENCE_SHEET_TO_MANIFEST.out.manifests
        .flatMap { manifest_dir ->
            manifest_dir.listFiles()
                .findAll { it.isDirectory() }
                .collect { genus_dir ->
                    def genus = genus_dir.name
                    def meta = [ id: genus.replaceAll('\\s+', '_'), genus: genus ]
                    [ meta, file("${genus_dir}/references.tsv"), file("${genus_dir}/ref_list.txt"), file("${genus_dir}/staged_refs") ]
                }
        }

    ch_reference_tables_from_ncbi = NCBI_DATASETS_DOWNLOAD.out.references
        .join(NCBI_DATASETS_DOWNLOAD.out.ref_list)
        .join(NCBI_DATASETS_DOWNLOAD.out.staged_refs)
        .map { meta, references, ref_list, staged_refs -> [ meta, references, ref_list, staged_refs ] }

    ch_reference_tables = ch_reference_tables_from_sheet
        .mix(ch_reference_tables_from_ncbi)

    //
    // FastANI: query vs type-strain references
    //
    ch_fastani_input = PREPARE_QUERY_LIST.out.query_list
        .join(ch_reference_tables)
        .map { meta, query_list, queries, references, ref_list, staged_refs ->
            [ meta, query_list, queries, ref_list, staged_refs ]
        }

    FASTANI_QUERY_VS_REFERENCE(
        ch_fastani_input,
    )
    ch_versions = ch_versions.mix(FASTANI_QUERY_VS_REFERENCE.out.versions)

    //
    // Optional intra-genus all-vs-all ANI clustering
    //
    ch_intra_input = params.skip_intra_genus_ani
        ? channel.empty()
        : PREPARE_QUERY_LIST.out.query_list

    FASTANI_INTRA_GENUS(
        ch_intra_input,
    )
    ch_versions = ch_versions.mix(FASTANI_INTRA_GENUS.out.versions)

    ch_intra_ani = params.skip_intra_genus_ani
        ? channel.empty()
        : FASTANI_INTRA_GENUS.out.ani

    //
    // Optional ANI heatmaps per genus
    //
    ch_heatmap_input = params.skip_heatmap
        ? channel.empty()
        : FASTANI_QUERY_VS_REFERENCE.out.ani.mix(ch_intra_ani)

    ANI_HEATMAP(
        ch_heatmap_input,
    )
    ch_versions = ch_versions.mix(ANI_HEATMAP.out.versions)

    //
    // Digital DDH on top-N references from FastANI
    //
    ch_all_references = ch_reference_tables
        .map { _meta, references, _ref_list, _staged_refs -> references }
        .collect()

    ch_select_ddh_ani = params.skip_ddh
        ? channel.empty()
        : FASTANI_QUERY_VS_REFERENCE.out.ani.map { _meta, tsv -> tsv }.collect()

    SELECT_DDH_PAIRS(
        ch_select_ddh_ani,
        ch_all_references,
        ASSIGN_GENUS.out.assignments,
    )
    ch_versions = ch_versions.mix(SELECT_DDH_PAIRS.out.versions)

    ch_ddh_pairs = params.skip_ddh
        ? channel.empty()
        : SELECT_DDH_PAIRS.out.pairs

    ch_empty_ddh_pairs = file("${projectDir}/assets/empty_ddh_pairs.tsv")

    ch_ddh_pairs_for_decision = params.skip_ddh
        ? channel.value(ch_empty_ddh_pairs)
        : ch_ddh_pairs

    ch_ddh_pairs_for_local = !params.skip_ddh && params.ddh_method in ['local', 'pyani']
        ? ch_ddh_pairs
        : channel.empty()

    ch_query_fasta_by_sample = PREPARE_QUERY_LIST.out.query_list
        .flatMap { _meta, _query_list, queries_dir ->
            queries_dir.listFiles()
                .findAll { it.name.endsWith('.fna') }
                .collect { fasta -> [ fasta.baseName, fasta ] }
        }

    ch_ref_fasta_by_accession = ch_reference_tables
        .flatMap { _meta, _references, _ref_list, staged_refs ->
            staged_refs.listFiles()
                .findAll { it.name.endsWith('.fna') }
                .collect { fasta -> [ fasta.baseName, fasta ] }
        }

    ch_ddh_jobs = ch_ddh_pairs_for_local
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def meta = [
                id: row.sample,
                genus: row.genus,
                ref_accession: row.ref_accession,
                ref_organism: row.ref_organism,
                single_end: true,
            ]
            [ row.sample, row.ref_accession, meta ]
        }
        .combine(ch_query_fasta_by_sample, by: 0)
        .map { sample, ref_accession, meta, query_fasta ->
            [ ref_accession, meta, query_fasta ]
        }
        .combine(ch_ref_fasta_by_accession, by: 0)
        .map { _ref_accession, meta, query_fasta, reference_fasta ->
            [ meta, query_fasta, reference_fasta ]
        }

    GBDP_DDH(
        ch_ddh_jobs,
    )
    ch_versions = ch_versions.mix(GBDP_DDH.out.versions)

    ch_ggdc_pairs = params.skip_ddh || params.ddh_method != 'ggdc'
        ? channel.empty()
        : ch_ddh_pairs

    GGDC_BUNDLE(
        ch_ggdc_pairs,
    )
    ch_versions = ch_versions.mix(GGDC_BUNDLE.out.versions)

    ch_ggdc_ingest_results = params.ggdc_results && !params.skip_ddh && params.ddh_method == 'ggdc'
        ? channel.value(file(params.ggdc_results, checkIfExists: true))
        : channel.empty()

    ch_ggdc_ingest_pairs = params.ggdc_results && !params.skip_ddh && params.ddh_method == 'ggdc'
        ? ch_ddh_pairs
        : channel.empty()

    GGDC_INGEST(
        ch_ggdc_ingest_results,
        ch_ggdc_ingest_pairs,
    )
    ch_versions = ch_versions.mix(GGDC_INGEST.out.versions)

    ch_ddh_collected = channel.empty()
    if (params.skip_ddh) {
        ch_ddh_collected = channel.value([])
    }
    else if (params.ddh_method in ['local', 'pyani']) {
        ch_ddh_collected = GBDP_DDH.out.ddh.map { _meta, tsv -> tsv }.collect()
    }
    else if (params.ddh_method == 'ggdc' && params.ggdc_results) {
        ch_ddh_collected = GGDC_INGEST.out.ddh.map { _meta, tsv -> tsv }.collect()
    }
    else {
        ch_ddh_collected = channel.value([])
    }

    SPECIES_DECISION(
        ASSIGN_GENUS.out.assignments,
        FASTANI_QUERY_VS_REFERENCE.out.ani.map { _meta, tsv -> tsv }.collect(),
        ch_ddh_collected,
        ch_ddh_pairs_for_decision,
    )
    ch_versions = ch_versions.mix(SPECIES_DECISION.out.versions)

    emit:
    summary     = SPECIES_DECISION.out.summary
    report      = SPECIES_DECISION.out.report
    assignments = ASSIGN_GENUS.out.assignments
    versions    = ch_versions
}
