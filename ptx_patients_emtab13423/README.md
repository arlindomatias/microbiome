# Gut microbiota of breast cancer patients with paclitaxel-induced neuropathy (E-MTAB-13423)

Snakemake + R re-analysis of 16S rRNA V3-V4 amplicon data from Zhong et al. (2023), *Blockade of CCR5 suppresses
paclitaxel-induced peripheral neuropathic pain caused by increased deoxycholic acid*, *Cell Reports* 42:113386
(ArrayExpress [E-MTAB-13423](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13423),
ENA [ERP152473](https://www.ebi.ac.uk/ena/browser/view/ERP152473)). The rat experiment of the same paper is
re-analysed with the same pipeline in [`../ptx_rats_emtab13426`](../ptx_rats_emtab13426).

**Report:** [Gut microbiota of breast cancer patients with paclitaxel-induced neuropathy](https://arlindomatias.github.io/microbiome/ptx_patients_emtab13423/reports/report.html)

The authors deposited only a relative-abundance table of the 20 most abundant genera. This pipeline starts from
the raw reads and produces ASV **counts**, which allow count-based methods (Chao1, DESeq2, ANCOM-BC2,
rarefaction).

## Design

Feces of breast cancer patients who completed paclitaxel treatment, graded for peripheral neuropathy (PIPN) with
NCI-CTCAE v4.0 and EORTC QLQ-CIPN20, and of healthy controls. One sample per person.

| Group | Subjects | Samples analysed | Sequencing run |
|---|---|---|---|
| Healthy control | C1-C20 | 20 | GP (NovaSeq 2×250) |
| Low PIPN (grade 1) | L1-L30 | 29 (L14 excluded) | GP (NovaSeq 2×250) |
| High PIPN (grade ≥ 2) | H1-H30 | 30 | GP (NovaSeq 2×250) |

**Caveats:**

- **Not a paclitaxel effect alone.** The controls have no cancer and received no chemotherapy, and breast cancer
  regimens usually combine paclitaxel with other drugs. Low vs high PIPN is the comparison closest to the question
  of the paper.
- **Cross-sectional design.** There is no sample before treatment, and diet, age, BMI and other drugs are not in
  ArrayExpress, so they cannot be adjusted for.
- **L14 is excluded.** Its deposited R2 file (`ZTPSN21GP644-L14_R2.fq.gz`) is an exact copy of R1 (same md5), so
  the pairs cannot be merged. The sample stays in the read tracking table.

## Workflow

```
ENA submitted FASTQ ──► cutadapt (both mate orientations) ──► DADA2 filterAndTrim
                    ──► DADA2 learnErrors (monotonic loess)/dada ──► mergePairs ──► removeBimeraDenovo
                    ──► SILVA 138.2 assignTaxonomy + addSpecies
                    ──► MAFFT + FastTree (midpoint-rooted) ──► phyloseq
                    ──► analysis: composition, alpha/beta diversity, PERMANOVA, DESeq2 + ANCOM-BC2 ──► report
```

The rules and scripts are the same as in the rat project (`../scripts/amplicon/`); the study design (groups,
comparisons, subject column, excluded samples, taxa of the paper) comes from [`config.yaml`](config.yaml).

| Rule | Script (in `../scripts/`) | What it does |
|---|---|---|
| `download_fastq`, `download_silva` | (shell) | Submitted FASTQ from ENA and SILVA 138.2 for DADA2, md5-checked, written under `$DATA` |
| `sample_metadata` | (Snakefile) | `results/samples.tsv` from the ArrayExpress SDRF |
| `cutadapt` | `amplicon/trim_primers.sh` | Primer removal in both mate orientations |
| `quality_profiles`, `filter_trim` | `amplicon/quality_profiles.R`, `amplicon/filter_trim.R` | DADA2 quality profiles and filtering |
| `dada2_denoise` | `amplicon/dada2_denoise.R` | Error models, ASVs, merging, chimera removal, read tracking |
| `assign_taxonomy` | `amplicon/assign_taxonomy.R` | SILVA 138.2 taxonomy (in chunks, see below) |
| `align_asvs`, `build_tree` | (MAFFT, FastTree) | Midpoint-rooted ASV tree for UniFrac |
| `build_phyloseq` | `amplicon/build_phyloseq.R` | phyloseq object and ASV count table, without excluded samples |
| `analysis` | `amplicon/analysis.R` | Statistics, result tables and figures |
| `report` | `amplicon/report.py` | Self-contained HTML report (`html_report.py`, shared with bulk-rna-seq and proteomics) |

Quirks of the deposited data, handled explicitly:

- **R1 and R2 are separate ENA runs.** The ENA-generated FASTQ files carry different run accessions in their read
  names, and some of them cannot be downloaded. The pipeline downloads the **submitted** files instead, which keep
  the original Illumina read names (same sequences).
- **Binned quality scores.** The read headers show an Illumina NovaSeq run (ENA lists MiSeq) with only four
  quality values (Q 2, 11, 25, 37). The DADA2 error model uses a loess fit forced to be non-increasing with
  quality (`error_function: loess_monotonic`).
- **Memory.** 15,754 ASVs do not fit a single `assignTaxonomy` call in 15 GB of RAM. The ASVs are classified in
  chunks of 2,000 (4 threads), each in a forked process and saved as a checkpoint under `$DATA`
  (`taxonomy_chunk_size`, `taxonomy_threads`).

### Methods

- **ASVs:** reads truncated to 220/229 bp after primer removal (~20 bp overlap on the ~429 bp amplicon); both
  primers found in 92–97% of pairs. Consensus chimera removal (40% of merged reads). 14,776–99,429 reads per
  sample remain.
- **Filtering:** bacteria only (no Archaea, chloroplast or mitochondria). Rarefied to the smallest library
  (14,776 reads) for alpha diversity, Jaccard and unweighted UniFrac.
- **Group tests:** Wilcoxon rank-sum for each pair of groups, BH-adjusted.
- **Beta diversity:** Bray-Curtis, Jaccard, unweighted and weighted UniFrac; PERMANOVA (`adonis2`, 9,999
  permutations) and `betadisper`.
- **Differential abundance (genus):** DESeq2 (`poscounts`, `~ group`) and ANCOM-BC2, on genera present in ≥ 25%
  of the compared samples. Genera with q < 0.05 in both methods are treated as robust.

All parameters (paths, primers, truncation, prevalence, cutoffs, permutations, seed) are in
[`config.yaml`](config.yaml).

### Current results

- 79 samples, 15,711 bacterial ASVs, 857 genera; 98% of the reads have a genus assignment.
- **Patients vs healthy controls:** community composition differs (PERMANOVA Bray-Curtis R² = 0.08 for low and
  0.07 for high PIPN, padj < 0.001), with larger dispersion in patients (`betadisper` p < 0.001) and no change in
  alpha diversity. Controls have more Actinomycetota (20% vs 7% in low PIPN, padj = 0.007), mostly
  *Bifidobacterium* (18% vs 5%); the mean share of Pseudomonadota is higher in low PIPN (15% vs 2%, mostly
  *Escherichia-Shigella*), but driven by a few samples (Wilcoxon padj = 0.19).
- **Robust genera** (q < 0.05 in DESeq2 and ANCOM-BC2): 20 for control vs low PIPN (e.g. *Akkermansia*,
  *Segatella*, *Citrobacter*, *[Ruminococcus] gnavus* group, *Veillonella* up; *Megasphaera*, *Desulfovibrio*,
  *Bifidobacterium*, *Ruminococcus*, *Blautia* down) and 7 for control vs high PIPN (*Enterococcus*, *Rothia*,
  *Segatella*, *Limosilactobacillus* up; *Megasphaera*, *Desulfovibrio* down).
- **Low vs high PIPN:** no difference in alpha diversity, Bray-Curtis (R² = 0.02, padj = 0.10) or weighted
  UniFrac; a small presence/absence difference (Jaccard R² = 0.02, padj = 0.006). 5 robust genera, all higher in
  high PIPN (*Acinetobacter*, *Parasutterella*, *Lactococcus*, *Pseudoduganella*, *[Eubacterium] eligens* group).
- **Taxa of the paper:** *Bifidobacterium* is higher in controls (vs low PIPN, padj = 0.004), as reported.
  *Ruminococcus* and *Dialister* are lower in both patient groups than in controls. None of *Clostridium*,
  *Dorea*, *Ruminococcus* or *Dialister* differs between low and high PIPN; *Sutterella* is too rare (< 25% of
  samples) to be tested.
- Genus abundances agree with the authors' processed table (18 genera matched, Spearman rho = 0.79); no sign of
  swapped sample labels.

## Folders

| Folder | Contents |
|---|---|
| `data/` | SDRF and authors' genus table (ArrayExpress), ASV sequences and tree; `seqtab_nochim.rds`, `ps.rds` and the alignment are not versioned |
| `results/` | Result tables (`.tsv`): samples, read tracking, taxonomy, ASV counts, statistics |
| `plots/` | All figures (PNG, embedded in the report) |
| `envs/` | Conda environments (`<name>.yaml` creates the environment `<name>`) and their `<name>.linux-64.pin.txt` locks |
| `reports/` | `report.html` (self-contained) and `sessionInfo.txt` |

FASTQ files, primer-trimmed and filtered reads, taxonomy checkpoints (`work/`) and SILVA are stored under `$DATA`
(see `config.yaml`).

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
