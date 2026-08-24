# nf-core/novelspecies: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/novelspecies, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

- NCBI type-strain download now uses `--from-type` and `--limit`, and resolves the *Bacillus* name collision (bacteria taxid 1386 vs walking-stick insects taxid 55087). GTDB suffixes (`Bacillus_A`, `Bacillus_AE`) are stripped only for the NCBI lookup.

### `Dependencies`

### `Deprecated`
