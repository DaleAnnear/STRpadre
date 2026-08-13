<p align="center">
  <img src="docs/cef4f14c-eb0e-4b92-8a6b-6c4cb90c1caa.png" width="260" alt="STRpadre logo">
</p>

<h1 align="center">STRpadre</h1>

<p align="center">
  A containerised Nextflow workflow for targeted tandem-repeat genotyping from PacBio HiFi and Oxford Nanopore long-read alignments.
</p>

<div align="center">

[![CI](https://github.com/DaleAnnear/STRpadre/actions/workflows/ci.yml/badge.svg)](https://github.com/DaleAnnear/STRpadre/actions/workflows/ci.yml)
[![GitHub issues](https://img.shields.io/github/issues/DaleAnnear/STRpadre.svg)](https://github.com/DaleAnnear/STRpadre/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/DaleAnnear/STRpadre.svg)](https://github.com/DaleAnnear/STRpadre/pulls)
[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-daleannear%2Fstrpadre-2496ED?logo=docker)](https://hub.docker.com/r/daleannear/strpadre)
[![Status](https://img.shields.io/badge/status-research%20use%20only-orange.svg)]()

</div>

---

<p align="center">
  Run complementary repeat genotypers, preserve their native output, and produce a traceable consensus for a predefined locus set.
</p>

> **Research use only  -  not for diagnostic use.** Tandem-repeat calls, caller
> agreement, and the STRpadre consensus require locus- and laboratory-specific
> validation. They must not be used as a clinical result without an appropriate
> validated process.

## Table of contents

- [About](#about)
- [Workflow overview](#workflow-overview)
- [Supported genotypers](#supported-genotypers)
- [Getting started](#getting-started)
- [Inputs and configuration](#inputs-and-configuration)
- [Running STRpadre](#running-strpadre)
- [Outputs and consensus](#outputs-and-consensus)
- [Testing and troubleshooting](#testing-and-troubleshooting)
- [References and licences](#references-and-licences)
- [Maintainer](#maintainer)

## About <a name="about"></a>

STRpadre is a Nextflow DSL2 workflow for targeted tandem-repeat genotyping from
pre-aligned PacBio HiFi and Oxford Nanopore (ONT) BAM/CRAM files. It starts after
alignment: it does not align raw reads, modify read groups, or infer a reference
build. The workflow validates inputs before callers begin, retains every native
caller result, converts calls through caller-specific loss-aware normalisers, and
builds a deterministic, auditable consensus.

The workflow is deliberately fail-fast. A caller failure stops the run by default;
it is not silently reclassified as a no-call. This prevents a failed tool from
biasing caller-agreement or consensus statistics.

## Workflow overview <a name="workflow-overview"></a>

~~~text
sample sheet + BAM/CRAM + reference + canonical loci
                         |
                         v
                containerised fail-fast preflight
                         |
   +---------------------+---------------------+---------------------+
   |                     |                     |                     |
 TRGT                  LongTR                ATaRVa               STRdust
 (HiFi)              (HiFi/ONT)            (HiFi/ONT)           (HiFi/ONT)
   |                     |                     |                     |
   +---- native VCFs, logs, and tool versions are retained ----------+
                         |
                         v
              caller-specific loss-aware normalisation
                         |
                         v
       deterministic allele assignment + consensus + QC reports
~~~

## Supported genotypers <a name="supported-genotypers"></a>

| Caller | Pinned version | Platforms | Default container | Native catalog | Key requirement |
|---|---:|---|---|---|---|
| [TRGT](https://github.com/PacificBiosciences/trgt) | 5.1.0 | HiFi | <code>local/strpadre-trgt:5.1.0</code> | structured BED (<code>ID;MOTIFS;STRUC</code>) | Locally built from a separately authorised PacBio binary |
| [LongTR](https://github.com/gymrek-lab/LongTR) | 1.2 | HiFi, ONT | <code>daleannear/strpadre:longtr-1.2-r1</code> | 1-based-start BED, motif, optional name | Coordinate-sorted/indexed BAM/CRAM with valid read-group metadata |
| [ATaRVa](https://github.com/SowpatiLab/ATaRVa) | 0.7.1 | HiFi, ONT | <code>daleannear/strpadre:atarva-0.7.1-r1</code> | sorted bgzip+tabix BED | Coordinate-sorted BAM/CRAM; MD/CS or <code>=/X</code> CIGAR is preferred |
| [STRdust](https://github.com/wdecoster/STRdust) | 0.20.0 | HiFi, ONT | <code>daleannear/strpadre:strdust-0.20.0-r1</code> | standard BED | BAM/CRAM plus reference; choose the unphased policy explicitly |

The caller configuration files in [configs/](configs) define the image, catalog,
platform support, resource request, and safe caller options. Exact source revisions
and build provenance are recorded in
[containers/versions.yml](containers/versions.yml).

## Getting started <a name="getting-started"></a>

### Prerequisites

Run Nextflow and Docker in the same Linux environment. On Windows, WSL2 with Docker
Desktop integration is the usual setup. You need:

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) **24.10 or later**
- [Docker Engine](https://docs.docker.com/engine/install/) or Docker Desktop with
  Linux/WSL integration enabled
- Git, Bash, and Python 3
- Coordinate-sorted and indexed BAM or CRAM input files
- A reference FASTA and matching <code>.fai</code>, one consistent reference build,
  a canonical locus manifest, and caller-native catalogs

Confirm the required runtimes before proceeding:

~~~bash
nextflow -version
docker --version
docker run --rm hello-world
~~~

### 1. Clone STRpadre

~~~bash
git clone https://github.com/DaleAnnear/STRpadre.git
cd STRpadre
~~~

### 2. Obtain the standard runtime images

LongTR, ATaRVa, STRdust, and the normaliser/consensus image are published on
[Docker Hub](https://hub.docker.com/r/daleannear/strpadre). Docker pulls them
automatically on the first run, or you can fetch them ahead of time:

~~~bash
docker pull daleannear/strpadre:normalizer-1.0.0-r1
docker pull daleannear/strpadre:longtr-1.2-r1
docker pull daleannear/strpadre:atarva-0.7.1-r1
docker pull daleannear/strpadre:strdust-0.20.0-r1
~~~

These are the default runtime images in <code>nextflow.config</code>. The caller YAML files retain image metadata for validation and provenance; use the matching <code>*_container</code> parameter to override a runtime image.

### 3. Build TRGT locally when you need it

TRGT is different from the other callers. The repository does **not** download,
bundle, or publish the TRGT binary. You must first obtain an authorised **TRGT
5.1.0** binary directly from PacBio under the
[PacBio Software License Agreement](https://github.com/PacificBiosciences/trgt/blob/main/LICENSE.md).

This local build is intentional:

- The PacBio agreement applies to the TRGT software; STRpadre cannot accept it on
  your behalf or redistribute a binary supplied under your organisation's agreement.
- The workflow uses the binary only in the local
  <code>local/strpadre-trgt:5.1.0</code> image. It is not included in the Docker Hub
  images.
- TRGT is restricted to PacBio HiFi input in this workflow. The preflight rejects
  TRGT for ONT samples.
- Under the PacBio agreement, use TRGT only to process or analyse data generated
  on a PacBio instrument or otherwise provided by PacBio.
- Do not commit the binary or publish an image containing it unless your organisation
  has confirmed that its PacBio agreement permits that distribution.

After accepting the agreement and obtaining the correct binary, place it at the
exact required path and build:

~~~bash
mkdir -p containers/vendor
install -m 0555 /path/to/authorised/trgt-5.1.0 containers/vendor/trgt-5.1.0

# Builds the local TRGT image and refreshes the other local build images.
bash containers/build.sh

# Confirms that every locally built image starts and reports a version.
bash containers/smoke_test.sh
~~~

[containers/Dockerfile.trgt](containers/Dockerfile.trgt) copies this file into the
image and checks that <code>trgt --version</code> reports <code>5.1.0</code>; the
build therefore fails rather than producing an image with an unexpected binary. This
is a practical summary, not legal advice: read the agreement that accompanied your
binary and follow your organisation's licensing process.

## Inputs and configuration <a name="inputs-and-configuration"></a>

### Prepare a sample sheet

Start from the checked-in templates:

~~~bash
cp examples/samples.csv samples.csv
cp examples/params.yml params.yml
~~~

Replace the template paths with absolute paths to your own data. The sample sheet must
contain the following columns; <code>sex</code> and <code>ploidy</code> are optional:

~~~csv
sample_id,platform,alignment,alignment_index,sex,ploidy
HG002,hifi,/data/HG002.sorted.bam,/data/HG002.sorted.bam.bai,male,2
NA12878,ont,/data/NA12878.sorted.cram,/data/NA12878.sorted.cram.crai,female,2
~~~

- <code>sample_id</code> must be unique and use only <code>A-Za-z0-9_.-</code>.
- <code>platform</code> is exactly <code>hifi</code> or <code>ont</code>.
- <code>alignment_index</code> must be the matching BAM/CRAM index.
- CRAM input must be decodable with the FASTA specified in <code>params.yml</code>.
- When the manifest includes <code>chrX</code> or <code>chrY</code> and LongTR or STRdust is selected, <code>sex</code> is required and must be <code>male</code>/<code>XY</code> or <code>female</code>/<code>XX</code>. The default GRCh38 policy calls male non-PAR X/Y loci haploid, male PAR loci diploid, female X loci diploid, and female Y loci as no-calls.
  For another assembly or contig naming scheme, update the <code>sex_chromosome_policy</code> in both caller configuration files to matching PAR intervals before running.

### Prepare the locus manifest and caller catalogs

The canonical locus manifest is TSV and uses **0-based, half-open** coordinates:

~~~tsv
locus_id	chr	start	end	coordinate_system	reference_build	primary_motif	repeat_structure	trgt_id	longtr_id	atarva_id	strdust_id
HTT	chr4	3074876	3074966	0-based-half-open	GRCh38	CAG		HTT	HTT	HTT	chr4:3074877-3074966
~~~

Read [docs/locus-manifest.md](docs/locus-manifest.md) before creating a manifest.
It is an identity and mapping layer, not a substitute for a caller's native catalog.

For simple loci, STRpadre can generate LongTR, ATaRVa, and STRdust catalogs and a
limited simple-repeat TRGT catalog:

~~~bash
nextflow run catalogs.nf -profile docker \
  --locus_manifest loci.tsv \
  --reference_fai reference.fa.fai \
  --callers longtr,atarva,strdust \
  --outdir catalog-build
~~~

Catalogues are generated from the supplied locus manifest and reference FAI. See
[docs/catalog-adapter.md](docs/catalog-adapter.md) for their caller-specific formats.
ATaRVa includes its corresponding <code>.bed.gz.tbi</code> index automatically.

### Configure <code>params.yml</code>

At minimum, set these fields to paths for your analysis:

~~~yaml
samplesheet: /absolute/path/samples.csv
reference: /absolute/path/reference.fa
reference_fai: /absolute/path/reference.fa.fai
locus_manifest: /absolute/path/loci.tsv
outdir: results
callers: [longtr, atarva, strdust]
~~~

Caller catalogues are generated automatically from `--locus-manifest` and `--reference-fai` for every supported genotyper. Caller YAML files configure the image and caller options only; `reference_build` is optional caller metadata. The canonical manifest remains the source of reference-build provenance.

<code>additional_args</code> is an argv list, not a shell fragment. It rejects
whitespace, shell metacharacters, and caller input/output flags owned by the workflow.

## Running STRpadre <a name="running-strpadre"></a>

### Run a lightweight workflow test

This checks Nextflow wiring without executing the native caller tools:

~~~bash
nextflow run main.nf -profile test -stub-run \
  -params-file tests/data/stub.params.yml
~~~

### Run a standard non-TRGT analysis

This is the usual command for an ONT cohort, or for a HiFi cohort that does not use
TRGT:

~~~bash
nextflow run main.nf -profile docker \
  -params-file params.yml \
  --callers longtr,atarva,strdust \
  -resume
~~~

### Run a HiFi analysis with TRGT

Only run this after completing the local TRGT build above:

~~~bash
nextflow run main.nf -profile docker \
  -params-file params.yml \
  --callers trgt,longtr,atarva,strdust \
  -resume
~~~

For a mixed HiFi/ONT cohort, the same caller list is valid: TRGT runs on the HiFi
rows only, and the caller matrix records TRGT as not applicable for ONT rows.

<code>-resume</code> reuses work only when inputs, code, configuration, and container
identity match. Resource defaults are process labels in
[nextflow.config](nextflow.config); use an institutional Nextflow profile to adjust
them. Do not change a caller's <code>threads</code> setting without consulting that
caller's documentation - LongTR is deliberately fixed at one calling thread.

LongTR uses <code>--lib-from-samp</code> by default, so each sample BAM/CRAM is treated
as one library and <code>@RG</code> <code>LB</code> tags are not required. Set
<code>options.use_lb_tags: true</code> in the LongTR YAML configuration only when every
read group carries a valid <code>LB</code> tag and library-aware calling is desired.

## Outputs and consensus <a name="outputs-and-consensus"></a>

With <code>--outdir results</code>, STRpadre creates predictable,
collision-free directories:

~~~text
results/
  validation/              preflight report, locus_mapping.tsv, versions.yml
  raw/<caller>/<sample>.{caller}.native/
                            native caller output, sorted/indexed calls.vcf.gz,
                            command.json, stderr log, versions.yml
  normalized/<caller>/      common loss-aware per-sample TSV.gz and parser versions
  consensus/
    consensus.tsv.gz        authoritative consensus per sample and locus
    consensus.vcf.gz(.tbi)  standards-compliant companion VCF
    caller_matrix.tsv.gz    selected/applicable/missing/filtered/no-call matrix
    discordant_calls.tsv.gz calls with discordance or inadequate support
    locus_summary.tsv.gz    cohort status counts per canonical locus
    provenance.json         input hashes and algorithm provenance
    qc_report.md            concise cohort QC summary
  reports/                  Nextflow report, trace, timeline, and DAG
~~~

The authoritative consensus is <code>consensus.tsv.gz</code>. The companion VCF
deliberately uses <code>GT=./.</code> and custom <code>AL</code>/<code>CN</code>
values rather than inventing one universal ALT-indexed genotype for tools that use
different sequence encodings.

Calls are joined through the canonical locus map, never by row order. For diploid
calls, STRpadre treats alleles as unordered unless reliable phase metadata is present.
The default consensus policy requires two supporting callers within 5 bp or 5%
relative length and does not allow filtered calls to support a consensus. It has no
sequence-level consensus. Review and archive
[configs/consensus.yml](configs/consensus.yml) with every analysis.

Important limitations:

- Copy number is retained only when a caller provides it or a full repeat length
  divides exactly by the canonical primary motif.
- Complex and interrupted repeats cannot safely be reduced to one copy-number value.
- TRGT padding and caller-specific sequence definitions are not assumed comparable.
- For exploratory partial analyses, rerun the successfully completed samples as a new,
  explicitly declared analysis rather than converting a caller failure into no-calls.

## Testing and troubleshooting <a name="testing-and-troubleshooting"></a>

Run the local code checks:

~~~bash
pytest -q
ruff check bin tests
python -m py_compile bin/*.py
~~~

See [docs/testing.md](docs/testing.md) for Docker integration and Nextflow stub
tests. Common causes of a failed preflight are:

- a reference-build or contig-name mismatch;
- missing reference access for a CRAM file;
- unsorted or unindexed alignments;
- missing LongTR <code>@RG</code> sample/library fields; or
- an ATaRVa BED that is not bgzip-compressed and tabix-indexed.

## References and licences <a name="references-and-licences"></a>

Please cite the caller publications that support your analysis:

1. **TRGT**  -  Dolzhenko E *et al.* Characterization and visualization of tandem
   repeats at genome scale. *Nature Biotechnology* **42**, 1606-1614 (2024).
   [doi:10.1038/s41587-023-02057-3](https://doi.org/10.1038/s41587-023-02057-3)
2. **LongTR**  -  Jam H, Zook J, Javadzadeh S, Park J, Sehgal A & Gymrek M.
   Genome-wide profiling of genetic variation at tandem repeats from long reads.
   *Genome Biology* **25**, 176 (2024).
   [doi:10.1186/s13059-024-03319-2](https://doi.org/10.1186/s13059-024-03319-2)
3. **ATaRVa**  -  Sivakumar AK *et al.* ATaRVa: Analysis of Tandem Repeat Variation
   from Long Read Sequencing data. *bioRxiv* (preprint, 2025).
   [doi:10.1101/2025.05.13.653434](https://doi.org/10.1101/2025.05.13.653434)
4. **STRdust**  -  De Coster W *et al.* Visualization and analysis of medically
   relevant tandem repeats in nanopore sequencing of control cohorts with pathSTR.
   *Genome Research* **34**, 2074-2080 (2024).
   [doi:10.1101/gr.279265.124](https://doi.org/10.1101/gr.279265.124)

LongTR is GPL-2.0-only; ATaRVa and STRdust are MIT-licensed. TRGT is governed by the
[PacBio Software License Agreement](https://github.com/PacificBiosciences/trgt/blob/main/LICENSE.md);
its separate binary and data-use restrictions are the reason for the local build
described above. Versions, upstream revisions, and container build facts are recorded
in [containers/versions.yml](containers/versions.yml).

## Maintainer <a name="maintainer"></a>

Dale J. Annear  -  [@DaleAnnear](https://github.com/DaleAnnear)
