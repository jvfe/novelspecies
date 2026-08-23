# nf-core/novelspecies: Genome-based species delimitation

[![Nextflow](https://img.shields.io/badge/nextflow-≥25.10.4-brightgreen.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/nf--core-1.0.0dev-brightgreen.svg)](https://nf-co.re/)

## Introduction

**nf-core/novelspecies** is a modular Nextflow DSL2 pipeline for determining whether bacterial isolate genomes or MAGs represent novel species relative to NCBI type strains.

The workflow:

1. Assigns a target genus to each query genome (GTDB-Tk summary, manual map, or samplesheet column)
2. Downloads genus-matched RefSeq type-strain assemblies from NCBI Datasets (or accepts a pre-built reference sheet)
3. Computes FastANI between queries and references, plus optional intra-genus all-vs-all ANI
4. Calculates digital DDH (Formula 2 / d4) for the top-*N* references per query
5. Applies standard species boundaries (ANI ≥ 95%, dDDH ≥ 70%) and generates summary tables plus HTML and heatmaps

## Quick start

```bash
nextflow run nf-core/novelspecies \
  -profile docker \
  --input samplesheet.csv \
  --gtdbtk_summary gtdbtk.bac120.summary.tsv \
  --outdir results/
```

### Samplesheet (`--input`)

```csv
sample,fasta
MT1012,/path/to/MT1012.fasta
MAG_bin23,/path/to/MAG_bin23.fna.gz
```

Optional `genus` column overrides GTDB-Tk/manual mapping for that sample.

### Taxonomy input (one required unless `genus` is in the samplesheet)

| Parameter | Description |
|-----------|-------------|
| `--gtdbtk_summary` | GTDB-Tk `bac120`/`ar53` summary TSV |
| `--genus_map` | CSV with `sample,genus` |

### Digital DDH modes

| `--ddh_method` | Behaviour |
|----------------|-----------|
| `local` (default) | Offline BLASTn GBDP formula-2 estimate |
| `pyani` | Alias for `local` (BLAST-based d4; not the pyani package) |
| `ggdc` | Builds manual GGDC web submission bundles; ingest results with `--ggdc_results` |

> **Note:** The official [GGDC web service](https://ggdc.dsmz.de) does not expose a stable public REST API. The `ggdc` mode prepares upload bundles and parses emailed/downloaded result files—it does not automate web submission.

## Outputs

| Path | Description |
|------|-------------|
| `species_delimitation_summary.tsv` | Per-sample verdict table |
| `summary/species_delimitation_report.html` | HTML dashboard |
| `ani/` | FastANI tables |
| `ani/heatmaps/` | Clustered ANI heatmaps per genus |
| `references/<genus>/` | Downloaded type-strain manifest |
| `ddh/pairwise/` | Local dDDH results (when enabled) |

### Verdict categories

- `same_species` — ANI ≥ 95% **and** dDDH ≥ 70%
- `novel_species_candidate` — ANI < 95% **and** dDDH < 70%
- `borderline_manual_review` — 94% ≤ ANI < 96%
- `ambiguous_mixed_signals` — conflicting ANI/dDDH signals

## Pre-computed references

Provide a reference sheet to skip NCBI download:

```csv
genus,accession,fasta,organism_name,is_type_strain
Klebsiella,GCF_000742135.1,/path/to/type.fna,Klebsiella pneumoniae,true
```

```bash
nextflow run nf-core/novelspecies \
  -profile apptainer \
  --input samplesheet.csv \
  --reference_sheet references.csv \
  --outdir results/
```

## Testing

Stub test (minimal synthetic genomes, dDDH skipped):

```bash
nextflow run . -profile test,docker
nf-test test tests/default.nf.test --profile test,docker
```

Integration test with real *Buchnera aphidicola* RefSeq genomes (~640 kb each), local dDDH, and species-boundary assertions:

```bash
nextflow run . -profile test_integration,docker
nf-test test tests/integration.nf.test --profile test_integration,docker
```

Test data live under `tests/data/integration/`; see `tests/data/integration/README.md` for accessions and expected outcomes.

## Credits

Built with [nf-core/tools](https://nf-co.re/tools) template v4.1.0.

## Citations

See [CITATIONS.md](CITATIONS.md).
