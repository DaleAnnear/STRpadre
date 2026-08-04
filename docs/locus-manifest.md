# Canonical locus manifest and native catalog adapters

The manifest is a tab-separated identity layer with one row per biological locus. It
uses 0-based, half-open intervals (`start` included, `end` excluded). Required columns
are `locus_id`, `chr`, `start`, `end`, `coordinate_system`, `reference_build`, and
`primary_motif`. Optional `repeat_structure` is descriptive only. The optional
`trgt_id`, `longtr_id`, `atarva_id`, and `strdust_id` columns explicitly map native
caller identifiers.

During preflight, every selected catalog record maps by native ID plus compatible
coordinates/motif, or by a unique coordinate/motif identity. It is never joined by
catalog row position. The resulting `validation/locus_mapping.tsv` is the consensus
join source.

Native formats are intentionally not interchangeable:

- TRGT requires a structured BED fourth column containing `ID`, `MOTIFS`, and `STRUC`.
  For complex loci, provide a native catalog.
- LongTR accepts a BED-like file whose start is documented as 1-based and whose fourth
  column is motif; its optional fifth column becomes the VCF ID.
- ATaRVa needs sorted, bgzip-compressed, tabix-indexed BED with motif and motif length.
- STRdust consumes a standard BED (zero-based half-open coordinates).

