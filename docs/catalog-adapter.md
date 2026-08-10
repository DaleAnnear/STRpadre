# Canonical catalog adaptation

`catalogs.nf` is the reusable Dockerized module for producing all selected caller
catalogs from **one canonical TSV locus manifest**. It is deliberately separate from
genotyping so that the generated files and their provenance can be reviewed, committed
to an analysis bundle, and referenced from the caller YAML files before a costly run.

```bash
nextflow run catalogs.nf -profile docker \
  --locus_manifest loci.tsv \
  --reference_fai reference.fa.fai \
  --callers trgt,longtr,atarva,strdust \
  --outdir catalog-build
```

Use `--callers all` for the same four-target selection. Names are normalized to lower
case; unknown, duplicated, and empty selections fail before the container launches.

The module creates `catalog-build/catalogs/` containing only the requested formats:

| Target | Generated file | Conversion |
|---|---|---|
| TRGT | `trgt.bed` | 0-based BED with `ID`, one `MOTIFS` value, and `STRUC=<TR>` |
| LongTR | `longtr.bed` | documented 1-based-start region BED with motif and locus ID |
| ATaRVa | `atarva.bed.gz` and `.tbi` | sorted 0-based BED with motif and motif length, bgzip/tabix indexed |
| STRdust | `strdust.bed` | 0-based three-column BED |

It also writes `catalog_locus_mapping.tsv` and `catalog_adaptation_report.json`.
These identify every generated native locus and record the target coordinate convention;
they are not substituted for the workflow's later preflight mapping report.

## Safety policy

The adapter validates manifest columns, a single non-empty reference build, DNA motifs,
unique native identifiers, reference-contig names, and interval bounds against the
provided `.fai`. Output is sorted in reference-index order and created atomically: a
failed conversion does not leave a partial catalog directory.

A canonical primary motif is sufficient for a simple repeat, but cannot encode a
complex TRGT motif set or structure. Therefore the module generates TRGT only when
`repeat_structure` is empty. If it is non-empty, generation stops with an actionable
error; provide a reviewed native TRGT catalog instead. This prevents flattening a
complex locus into a scientifically misleading simple repeat.

## Using generated files in genotyping

`main.nf` invokes this adapter automatically, then passes its catalogues to input validation and to each caller. The generated files are the sole caller-catalogue source; static `native_catalog` entries in caller YAML files are not used. Run `catalogs.nf` or `bin/catalog_adapter.py` directly only when you want to inspect or archive the generated catalogues.

`bin/catalog_adapter.py` is the underlying command-line program used by the Nextflow
module. It accepts the same manifest, FAI, caller selection, and output directory
arguments, but production use should prefer `catalogs.nf -profile docker` so `pysam`,
bgzip, and tabix behavior are pinned.
