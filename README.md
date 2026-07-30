# STRpadre

STRpadre is a research-use-only Nextflow DSL2 workflow for targeted tandem-repeat
genotyping from pre-aligned PacBio HiFi and Oxford Nanopore (ONT) BAM/CRAM files.
It deliberately starts after alignment: it does not align raw reads, alter read groups,
or infer a reference build.

> **Not for diagnostic use.** Tandem-repeat calls, caller agreement, and this
> consensus are research results that require locus- and laboratory-specific validation.
> In particular, TRGT is governed by the PacBio Software License Agreement and may
> only process PacBio-provided data; this workflow hard-rejects TRGT for ONT samples.

## What it does

```text
sample sheet + BAM/CRAM + reference + canonical loci
                         │
                         ▼
                containerized fail-fast preflight
                         │
       ┌─────────────────┼──────────────────┬─────────────────┐
       ▼                 ▼                  ▼                 ▼
     TRGT             LongTR             ATaRVa           STRdust
   (HiFi)          (HiFi/ONT)          (HiFi/ONT)       (HiFi/ONT)
       │                 │                  │                 │
       └── native VCFs, logs, executable versions are retained ┘
                         │
                         ▼
              caller-specific loss-aware normalizers
                         │
                         ▼
       deterministic allele assignment + robust consensus + QC
```

Every caller process, normalizer, preflight, and consensus process has a Docker
container. Processes fail immediately by default; partial consensus is deliberately
disabled (`allow_partial_results: false`) because treating a crashed caller as a
no-call would bias agreement statistics.

## Caller compatibility and contracts

| Caller | Pinned version | HiFi | ONT | Native catalog | Important contract |
|---|---:|:---:|:---:|---|---|
| TRGT | 5.1.0 | yes | no | structured BED (`ID;MOTIFS;STRUC`) | aligned HiFi BAM; PacBio license restriction |
| LongTR | 1.2 | yes | yes | 1-based-start BED, motif, optional name | coordinate-sorted/indexed BAM/CRAM and valid RG metadata |
| ATaRVa | 0.7.1 | yes | yes | sorted bgzip+tabix BED: chrom/start/end/motif/motif length | coordinate-sorted BAM/CRAM; MD/CS or `=/X` CIGAR preferred |
| STRdust | 0.20.0 | yes | yes | standard BED | BAM/CRAM plus reference; configure `--unphased` policy explicitly |

The contracts above are implemented from the upstream [TRGT CLI and VCF
documentation](https://github.com/PacificBiosciences/trgt/tree/main/docs),
[LongTR documentation](https://github.com/gymrek-lab/LongTR),
[ATaRVa README/CLI](https://github.com/SowpatiLab/ATaRVa), and
[STRdust README/CLI](https://github.com/wdecoster/STRdust). Exact upstream source
revisions, licenses, container tags, and build policy are in
[containers/versions.yml](containers/versions.yml).

## Prerequisites

- Nextflow 24.10 or later
- Docker Engine
- Coordinate-sorted, indexed BAM or CRAM. CRAM must be decodable with the supplied
  reference FASTA.
- A single known reference build represented consistently by FASTA, `.fai`, catalogs,
  canonical manifest, and alignments.
- Caller images built or pulled before a run. The workflow never builds an image.

Build images once:

```bash
bash containers/build.sh
bash containers/smoke_test.sh
```

`containers/Dockerfile.trgt` intentionally requires a locally supplied,
license-authorized `containers/vendor/trgt-5.1.0` binary. Do not publish an image
containing TRGT without confirming the PacBio agreement. The other Dockerfiles build
from the pinned upstream revisions in `containers/versions.yml`.

## Inputs

The CSV sample sheet must contain these columns. `sex` and `ploidy` are optional.

```csv
sample_id,platform,alignment,alignment_index,sex,ploidy
HG002,hifi,/data/HG002.sorted.bam,/data/HG002.sorted.bam.bai,male,2
NA12878,ont,/data/NA12878.sorted.cram,/data/NA12878.sorted.cram.crai,female,2
```

- `sample_id` is unique and restricted to `A-Za-z0-9_.-`; it becomes part of output
  paths and is checked before callers launch.
- `platform` is exactly `hifi` or `ont`.
- `alignment_index` must be the matching BAM/CRAM index.
- The preflight process checks sample-sheet shape, file existence, reference `.fai`,
  canonical intervals, caller schemas/catalogs, catalog-to-locus mappings, alignment
  sort order, and BAM/CRAM/reference contig compatibility.

The canonical locus manifest is tab-separated and uses **0-based, half-open**
coordinates exclusively in this release:

```tsv
locus_id	contig	start	end	coordinate_system	reference_build	primary_motif	repeat_structure	trgt_id	longtr_id	atarva_id	strdust_id
HTT	chr4	3074876	3074966	0-based-half-open	GRCh38	CAG		HTT	HTT	HTT	chr4:3074877-3074966
```

Read [docs/locus-manifest.md](docs/locus-manifest.md) before preparing catalogs. A
canonical manifest is an identity/mapping layer, not a substitute for every caller's
native catalog fields.

## Configuration and execution

Separate strict schemas validate global parameters and each YAML file:

- [nextflow_schema.json](nextflow_schema.json)
- [configs/trgt.yml](configs/trgt.yml), [configs/longtr.yml](configs/longtr.yml),
  [configs/atarva.yml](configs/atarva.yml), [configs/strdust.yml](configs/strdust.yml)
- [configs/consensus.yml](configs/consensus.yml)

`additional_args` is an argv list. It rejects whitespace, shell metacharacters, and
caller input/output flags owned by the workflow. It is never evaluated by a shell.

The checked-in [examples/params.yml](examples/params.yml) and
[examples/samples.csv](examples/samples.csv) are a complete template; replace the
absolute alignment paths and `example-build` catalogs/reference with your data.
ATaRVa requires the configured `.bed.gz.tbi`; generate it reproducibly from the
canonical manifest when appropriate:

```bash
python bin/build_catalogs.py --caller atarva --locus-manifest loci.tsv --output atarva.bed.gz
```

Use native catalogs whenever a simple conversion would lose structure. `build_catalogs.py`
only converts LongTR, ATaRVa, and STRdust formats directly; TRGT conversion is limited
to simple, single-motif loci and requires `--allow-simple-trgt`.

Run with a YAML caller list:

```bash
nextflow run main.nf -profile docker -params-file params.yml -resume
```

or override it conveniently on the command line:

```bash
nextflow run main.nf -profile docker -params-file params.yml --callers trgt,longtr,atarva,strdust -resume
```

Three examples (after replacing the template paths):

```bash
# 1. HiFi: all supported callers
nextflow run main.nf -profile docker -params-file examples/params.yml --callers trgt,longtr,atarva,strdust -resume

# 2. ONT: TRGT is intentionally excluded
nextflow run main.nf -profile docker -params-file ont.params.yml --callers longtr,atarva,strdust -resume

# 3. Mixed cohort: TRGT runs only on HiFi rows; the caller matrix records not_applicable_to_platform for ONT
nextflow run main.nf -profile docker -params-file mixed.params.yml --callers trgt,longtr,atarva,strdust -resume
```

Resource defaults are process labels in [nextflow.config](nextflow.config). Override
them with an institutional profile; do not change a caller's `threads` field without
checking its native documentation (LongTR is deliberately fixed at one calling thread).
`-resume` reuses work directories only when inputs, process code, configuration, and
container identity match.

## Outputs

`--outdir results` produces predictable, collision-free directories:

```text
results/
  validation/              preflight report, locus_mapping.tsv, versions.yml
  raw/<caller>/<sample>.{caller}.native/
                            complete native caller output, sorted/indexed calls.vcf.gz,
                            command.json, stderr log, versions.yml
  normalized/<caller>/      common loss-aware per-sample TSV.gz and parser versions
  consensus/
    consensus.tsv.gz        consensus sizes/copy numbers/statuses/provenance per sample+locus
    consensus.vcf.gz(.tbi)  standards-compliant companion VCF
    caller_matrix.tsv.gz    selected/applicable/missing/filtered/no-call matrix
    discordant_calls.tsv.gz calls with discordance or inadequate support
    locus_summary.tsv.gz    cohort status counts per canonical locus
    provenance.json         input hashes, policy, algorithm provenance
    qc_report.md            concise cohort QC summary
  reports/                  Nextflow report, trace, timeline, DAG
```

The VCF deliberately uses `GT=./.` and custom `AL`/`CN` values: a universal ALT-indexed
genotype is not scientifically defensible when callers encode sequences differently.
The authoritative consensus genotype representation is `consensus.tsv.gz`; the VCF
contains documented compatibility/support/status fields without inventing a sequence
allele.

## Consensus method and limits

Calls join on the canonical locus map, never row order. The implementation verifies
coordinates, reference build, and motifs (with configurable exact/rotation/reverse-
complement equivalence). Diploid pairs are treated as unordered unless a caller gives
trustworthy phase metadata. Each new pair is assigned to current allele clusters by the
minimum total absolute length difference; exact ties follow caller name then input path.
Per-cluster medians yield sizes. Support, min/max range, contributing callers,
concordance, and every exclusion reason are retained.

The default policy requires two supporting callers within 5 bp or 5% relative length,
does not let filtered calls contribute, and has no sequence-level consensus. Change only
the explicit fields in `configs/consensus.yml`, then archive that file with the run.

Known limitations:

- Copy number is retained only where a caller provides a documented value or where a
  full base-pair repeat length divides exactly by the canonical primary motif.
- Complex repeat structures and interrupted motifs cannot safely be reduced to a
  single copy-number value.
- Sequence consensus is intentionally disabled; TRGT padding and caller-specific
  sequence definitions are not assumed comparable.
- Caller failure is fatal. For an exploratory partial cohort, rerun successful samples
  as a new declared analysis rather than silently converting failed work into no-calls.

## Development and troubleshooting

```bash
pytest -q
ruff check bin tests
python -m py_compile bin/*.py
nextflow run main.nf -profile test -stub-run -params-file tests/data/stub.params.yml
```

See [docs/testing.md](docs/testing.md) for the Docker integration and Nextflow stub
tests. Common failures are reference-build/contig mismatches, missing CRAM reference
access, unsorted or unindexed inputs, missing LongTR `@RG` sample/library fields, and
an ATaRVa BED that is not bgzip+tabix indexed.

## Citations and licenses

Please cite the caller publications described by their upstream projects. LongTR is
GPL-2.0-only; ATaRVa and STRdust are MIT; TRGT has the PacBio Software License
Agreement. The authoritative source URLs, versions, revisions, and build facts are
machine-readable in [containers/versions.yml](containers/versions.yml).
