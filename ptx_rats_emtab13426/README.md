# Gut microbiota of paclitaxel-treated rats (E-MTAB-13426)

Snakemake + R re-analysis of 16S rRNA V3-V4 amplicon data from Zhong et al. (2023), *Blockade of CCR5 suppresses
paclitaxel-induced peripheral neuropathic pain caused by increased deoxycholic acid*, *Cell Reports* 42:113386
(ArrayExpress [E-MTAB-13426](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13426),
ENA [PRJEB67467](https://www.ebi.ac.uk/ena/browser/view/PRJEB67467)).

**Report:** [Gut microbiota of paclitaxel-treated rats](https://arlindomatias.github.io/microbiome/ptx_rats_emtab13426/reports/report.html)

The authors deposited only relative-abundance tables. This pipeline starts from the raw reads and produces ASV
**counts**, which allow count-based methods (Chao1, DESeq2, ANCOM-BC2, rarefaction).

## Design

| Group | Rats | Samples | Sequencing run |
|---|---|---|---|
| Control (vehicle) | C1-C10 | 10 | AG (MiSeq 2×300) |
| PTX week 2 (day 14) | P1-P8 | 8 | AG (MiSeq 2×300) |
| PTX week 4 (day 28) | P1-P8 (same rats) | 8 | BD (MiSeq 2×250) |

**Caveats:**

- All week-4 samples come from a separate sequencing run, so run effects are confounded with time. Only Control
  vs PTX week 2 is free of this confounding.
- The week 2 and week 4 labels of the raw data (used here) look swapped relative to the authors' processed table:
  in 15 of 16 PTX samples the best-matching authors' sample is the same rat at the other time point.

## Workflow

```
ENA FASTQ ──► cutadapt (both mate orientations) ──► DADA2 filterAndTrim
          ──► DADA2 learnErrors/dada per run × orientation batch ──► mergePairs ──► removeBimeraDenovo
          ──► SILVA 138.2 assignTaxonomy + addSpecies
          ──► MAFFT + FastTree (midpoint-rooted) ──► phyloseq
          ──► analysis: composition, alpha/beta diversity, PERMANOVA, DESeq2 + ANCOM-BC2 ──► report
```

| Rule | Script (in `../scripts/`) | What it does |
|---|---|---|
| `download_fastq`, `download_silva` | (shell) | ENA FASTQ and SILVA 138.2 for DADA2, md5-checked, written under `$DATA` |
| `sample_metadata` | (Snakefile) | `results/samples.tsv` from the ArrayExpress SDRF |
| `cutadapt` | `amplicon/trim_primers.sh` | Primer removal in both mate orientations |
| `quality_profiles`, `filter_trim` | `amplicon/quality_profiles.R`, `amplicon/filter_trim.R` | DADA2 quality profiles and filtering |
| `dada2_denoise` | `amplicon/dada2_denoise.R` | Error models, ASVs, merging, chimera removal, read tracking |
| `assign_taxonomy` | `amplicon/assign_taxonomy.R` | SILVA 138.2 taxonomy |
| `align_asvs`, `build_tree` | (MAFFT, FastTree) | Midpoint-rooted ASV tree for UniFrac |
| `build_phyloseq` | `amplicon/build_phyloseq.R` | phyloseq object and ASV count table |
| `analysis` | `amplicon/analysis.R` | Statistics, result tables and figures |
| `report` | `amplicon/report.py` | Self-contained HTML report (`html_report.py`, shared with bulk-rna-seq and proteomics) |

Two quirks of the deposited data are handled explicitly:

- **Mixed mate orientation:** in 14 of the 26 runs R1 starts with the reverse primer (806R).
  `trim_primers.sh` tries both orientations and re-orients the reads, keeping 97.5-99.3% of pairs per sample.
- **Two read lengths:** the 2×250 run gives ~225/231 bp reads after primer removal. Truncation is set to
  220/229 bp for all samples, which leaves ~20 bp of overlap on the ~429 bp amplicon.

### Methods

- **ASVs:** DADA2 error models per sequencing run × orientation batch; consensus chimera removal (about 55% of
  merged reads are chimeric, most likely from the library PCR). 6,528–15,249 reads per sample remain.
- **Filtering:** bacteria only (no Archaea, chloroplast or mitochondria). Rarefied to the smallest library for
  alpha diversity, Jaccard and unweighted UniFrac.
- **Group tests:** Wilcoxon rank-sum for Control vs each PTX time point, signed-rank (paired by rat) for week 2 vs
  week 4, BH-adjusted.
- **Beta diversity:** Bray-Curtis, Jaccard, unweighted and weighted UniFrac; PERMANOVA (`adonis2`, 9,999
  permutations, restricted within rats for the paired comparison) and `betadisper`.
- **Differential abundance (genus):** DESeq2 (`poscounts`; `~ group` or `~ rat + group`) and ANCOM-BC2, on genera
  present in ≥ 25% of the compared samples. Genera with q < 0.05 in both methods are treated as robust.

All parameters (paths, primers, truncation, prevalence, cutoffs, permutations, seed) are in
[`config.yaml`](config.yaml).

### Current results

- 26 samples, 2,594 bacterial ASVs, 192 genera; 71% of the reads have a genus assignment.
- **Control vs PTX week 2:** community composition differs (PERMANOVA Bray-Curtis R² = 0.22, padj < 0.001) without
  a change in alpha diversity. 4 genera are called by both DA methods, all higher in PTX: UCG-005,
  Christensenellaceae R-7 group, *Desulfovibrio* and *Ligilactobacillus*.
- **PTX week 4:** a larger shift (Bray-Curtis R² = 0.34 vs control), higher alpha diversity and less Bacillota
  (58% vs 76% in controls). 17 genera are called by both methods (e.g. *Alistipes*, *Bacteroides*, *Roseburia* up;
  *Lactobacillus*, HT002 down). This time point is confounded with the sequencing run.
- Genus abundances agree with the authors' processed table (40 genera matched, Spearman rho = 0.90).

## Folders

| Folder | Contents |
|---|---|
| `data/` | SDRF and authors' genus table (ArrayExpress), ASV sequences and tree; `seqtab_nochim.rds`, `ps.rds` and the alignment are not versioned |
| `results/` | Result tables (`.tsv`): samples, read tracking, taxonomy, ASV counts, statistics |
| `plots/` | All figures (PNG, embedded in the report) |
| `envs/` | Conda environments (`<name>.yaml` creates the environment `<name>`) and their `<name>.linux-64.pin.txt` locks |
| `reports/` | `report.html` (self-contained) and `sessionInfo.txt` |

FASTQ files, primer-trimmed and filtered reads (`work/`) and SILVA are stored under `$DATA` (see `config.yaml`).

## Running

Requirements: Linux (x86-64), Conda and Snakemake ≥ 9. Exact versions are in `envs/*.linux-64.pin.txt`
(R 4.5.3, DADA2 1.38.0, phyloseq 1.54.2, DESeq2 1.50.2, ANCOMBC 2.14.0, vegan 2.7-5, cutadapt 5.2).

```bash
export DATA=/path/to/large/storage    # FASTQ, work files and SILVA are written here
snakemake -n --use-conda              # dry run
snakemake --use-conda --cores 8 --resources downloads=4
```

The phyloseq object (`data/ps.rds`) and the ASV count table (`results/asv_counts.tsv`) can be reused for other
analyses.
