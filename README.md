# Infernal ncRNA Search Pipeline

A pipeline for discovering non-coding RNA (ncRNA) homologs in genomic assemblies using
[INFERNAL](http://eddylab.org/infernal/) covariance model searches, following the
approach of the **Geronimo** tool:

> Kilar A.M., Fajkus P., & Fajkus J. (2023). **GERONIMO: A tool for systematic retrieval of structural RNAs in a broad evolutionary context.**
> *GigaScience*, 12, giad080. https://doi.org/10.1093/gigascience/giad080

The repository contains scripts to download genome assemblies from NCBI, run `cmsearch`
across all genomes and models, and parse the results into annotated CSV tables enriched
with taxonomy and flanking sequence context.

---

## Features

- **Standalone, parameterized parser** — `parse_tables.py` runs independently on any
  INFERNAL `cmsearch` output via a simple CLI; no Snakemake or R required.
- **Strand-aware sequence extraction** — minus-strand hits are automatically
  reverse-complemented and returned in correct 5′→3′ orientation.
- **Configurable flanking context** — extract each hit plus ±N bp of surrounding
  genomic sequence (`-r`), with contig-boundary clamping.
- **Lightweight extraction** — pulls sequences directly from FASTA with `seqtk`,
  avoiding a separate BLAST-database build step.
- **Resilient NCBI access** — automatic retries for taxonomy lookups (100×) and
  genome downloads (3×) for robust large-batch runs.
- **One-step download + full taxonomy** — produces a single `taxonomy.csv` with the
  complete lineage (superkingdom → genus).
- **HPC-ready** — PBS Pro batch scripts with explicit resource requests, tested at
  genome-batch scale on MetaCentrum.

---

## Repository contents

```
.
├── 01_02_download_taxonomy.sh   # Download genomes + build taxonomy.csv from NCBI
├── 04_pip_infernal_model.sh     # Run cmsearch (adapted from GERONIMO — see Credits)
├── 05_run_parse_tables.sh       # PBS wrapper — set paths here, then submit
├── parse_tables.py              # Parse and integrate all cmsearch results
└── README.md
```

---

## Dependencies

| Tool | Notes |
|---|---|
| **INFERNAL** | provides `cmbuild`, `cmcalibrate`, `cmsearch` |
| **NCBI E-Direct** | provides `esearch`, `efetch`, `xtract` |
| **seqtk** | sequence extraction with flanking context |
| **Python 3.6+** | |
| **pandas** | `pip install pandas` |
| **wget**, **gunzip** | standard Unix utilities |

Scripts are written for a **PBS Pro** cluster environment (tested on MetaCentrum).
Typical resource requirements per job: 10–12 CPUs, 45 GB RAM, 20 h walltime.

---

## Pipeline

### Step 1 — Download genomes and build taxonomy table

```bash
bash 01_02_download_taxonomy.sh
```

Queries NCBI Assembly via E-Direct, retrieves the full taxonomic lineage
(superkingdom → genus) for each result, and downloads the genomic FASTA files.
Failed queries are retried up to 100 times; failed downloads up to 3 times.

Set the NCBI query string and output directory at the top of the script before running.

**Outputs:**
- `taxonomy.csv` — one row per assembly; columns: `superkingdom`, `kingdom`, `phylum`,
  `class`, `order`, `family`, `genus`, `ScientificName`, `GCA`, `taxid`
- `*.fna` — decompressed genomic FASTA files

---

### Step 2 — Run cmsearch

```bash
bash 04_pip_infernal_model.sh
```

Loops over every covariance model (`cov_*`) in `DIR_MODELS` and every `*.fna` genome
in `DIR_DATABASE`, running `cmsearch` for each pair.

Set `DIR_MODELS` and `DIR_DATABASE` at the top of the script. Covariance models must be
built and calibrated beforehand with `cmbuild` / `cmcalibrate` (not included here;
see the INFERNAL documentation or the Geronimo pipeline linked above).

**Outputs** (written to `DIR_DATABASE/<model_name>/` for each model):

| File | Contents |
|---|---|
| `result_MODEL_vs_GENOME` | Full pairwise cmsearch output |
| `result_MODEL_vs_GENOME-alignment` | Structure alignments (Stockholm) |
| `result_MODEL_vs_GENOME.csv` | Tabular hit summary |

---

### Step 3 — Parse and integrate results

Edit the path variables in `05_run_parse_tables.sh`, then submit or run it:

```bash
bash 05_run_parse_tables.sh
```

This calls `parse_tables.py` directly. You can also invoke the script manually:

```bash
python3 parse_tables.py \
  -i /path/to/cmsearch/results \
  -t /path/to/taxonomy.csv \
  -r 200 \
  -d /path/to/genome/fastas \
  -o /path/to/output/dir
```

| Argument | Description |
|---|---|
| `-i` | Directory containing cmsearch result files |
| `-t` | `taxonomy.csv` from step 1 |
| `-r` | Flanking region to extract around each hit (bp) |
| `-d` | Directory containing genomic FASTA files |
| `-o` | Output directory |

`parse_tables.py` performs:
1. Parses all tabular cmsearch output files
2. Separates high-confidence hits (`!`) from marginal hits (`?`)
3. Extracts alignment strings and secondary structure annotation from alignment files
4. Merges taxonomy metadata (joined on GCA accession)
5. Extracts ±`-r` bp of flanking genomic sequence for each hit using `seqtk`
6. Counts hits per genome × model pair

**Outputs:**

| File | Contents |
|---|---|
| `infernal_result_hits.csv` | High-confidence hits, sorted by E-value |
| `infernal_result_possible_hits.csv` | Marginal hits, sorted by E-value |
| `infernal_result_report.csv` | Hit counts per genome × model pair |

#### Output columns

| Column | Description |
|---|---|
| `target_name` | Sequence name in the genome FASTA |
| `query_name` | Covariance model name |
| `seq_from` / `seq_to` | Hit coordinates on the target sequence |
| `strand` | `+` or `-` |
| `score` | Bit score |
| `E-value` | INFERNAL E-value |
| `model` | RNA family (parsed from filename) |
| `GCA` | NCBI GenBank Assembly accession |
| `ScientificName` | Species name |
| `family` … `superkingdom` | Full taxonomic lineage |
| `alignment` | Hit alignment string from INFERNAL |
| `GC_SS_cons` | Consensus secondary structure annotation |
| `GC_RF` | Reference annotation from the covariance model |
| `seq_from_extend` / `seq_to_extend` | Extended (flanking) coordinates |
| `extend_sequence` | Genomic sequence including ±`-r` bp flanking context |
| `number` | Hit rank within each genome × model pair |
| `number_of_hits` / `number_of_possible_hits` | Counts per genome × model |
