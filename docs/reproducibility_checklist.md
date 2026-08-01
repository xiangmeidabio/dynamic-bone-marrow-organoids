# Pre-release reproducibility checklist

Complete this checklist before creating the GitHub release cited by the manuscript.

## Scientific metadata

- [ ] Add the final manuscript title, author list, and corresponding author.
- [ ] Add the journal citation, DOI, and preprint link when available.
- [ ] Add raw and processed data accessions to `data/README.md`.
- [ ] Verify the organism and genome build for every dataset.
- [ ] Cite all external reference datasets and record their licenses.
- [ ] Map every main and extended-data figure panel to a script/notebook and output file.

## Analysis verification

- [ ] Confirm all QC thresholds against the final Methods section.
- [ ] Confirm expected doublet rates against 10x loading metrics.
- [ ] Confirm the DecontX contamination threshold and document its rationale.
- [ ] Confirm PCA dimensions, clustering resolutions, and manual cluster mappings.
- [ ] Confirm that `org.Hs.eg.db` or `org.Mm.eg.db` matches the species used in each enrichment analysis.
- [ ] Confirm that all cell-type labels exactly match the manuscript.
- [ ] Run every script/notebook from a clean environment and an empty `data/processed/` directory.
- [ ] Compare all regenerated figures with the submitted manuscript figures.

## Environment and provenance

- [ ] Commit a validated `renv.lock`.
- [ ] Freeze exact Python dependencies.
- [ ] Record R, Python, Seurat, CellChat, pySCENIC, CellOracle, and scTenifoldKnk versions.
- [ ] Record pySCENIC database URLs, release identifiers, and checksums.
- [ ] Record SHA-256 checksums for deposited processed objects.
- [ ] Record hardware and peak-memory requirements for large integration steps.

## GitHub and archival release

- [ ] Choose a software license with the institution and all authors.
- [ ] Complete and rename `CITATION.cff.template` to `CITATION.cff`.
- [ ] Remove all `TODO` fields.
- [ ] Check that no patient identifiers, credentials, absolute paths, or local usernames remain.
- [ ] Create a versioned GitHub release matching the manuscript revision.
- [ ] Archive that release in Zenodo or an equivalent repository and add the concept DOI.
- [ ] Cite the immutable release DOI in the Code Availability statement.
