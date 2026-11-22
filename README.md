# Microbiome: chemotherapy-induced changes in the gut microbiota

Reproducible Snakemake + R workflows for 16S rRNA amplicon analysis of public microbiome datasets on
chemotherapy side effects. They go from the raw reads deposited in ENA to ASV counts, diversity, differential
abundance and an HTML report.

| Folder | Dataset | Input | What it does |
|---|---|---|---|
| [`ptx_rats_emtab13426/`](ptx_rats_emtab13426) | [E-MTAB-13426](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13426) | ENA FASTQ (PRJEB67467) | cutadapt → DADA2 → SILVA 138.2 → phyloseq → alpha/beta diversity, PERMANOVA, DESeq2 + ANCOM-BC2 |
| [`ptx_patients_emtab13423/`](ptx_patients_emtab13423) | [E-MTAB-13423](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13423) | ENA FASTQ (ERP152473) | Same pipeline: breast cancer patients with low or high paclitaxel-induced neuropathy vs healthy controls |

```
ENA FASTQ (PRJEB67467) ──► cutadapt ──► DADA2 ──► SILVA 138.2 + tree ──► phyloseq ──► statistics ──► report
```

The companion repositories [bulk-rna-seq](https://github.com/arlindomatias/bulk-rna-seq) and
[proteomics](https://github.com/arlindomatias/proteomics) follow the same conventions and report layout.
This repository was split out of [arlindomatias/projects](https://github.com/arlindomatias/projects), with its history.

## Reports

| Project | Report |
|---|---|
| `ptx_rats_emtab13426` | [Gut microbiota of paclitaxel-treated rats](https://arlindomatias.github.io/microbiome/ptx_rats_emtab13426/reports/report.html) |
| `ptx_patients_emtab13423` | [Gut microbiota of breast cancer patients with paclitaxel-induced neuropathy](https://arlindomatias.github.io/microbiome/ptx_patients_emtab13423/reports/report.html) |

Each report is built by the `report` rule of the project's Snakefile (`reports/report.html`) and published with GitHub Pages.

## Project layout

Every project follows the same folder layout. The scripts are not inside the projects: they live in the
repository-level [`scripts/`](scripts) folder, shared by the projects.

| Folder | Contents |
|---|---|
| `data/` | Inputs (sample sheet, authors' tables), ASV sequences and tree, intermediate objects (`seqtab_nochim.rds`, `ps.rds`, not versioned) |
| `results/` | Result tables only (`.tsv`) |
| `plots/` | All figures |
| `envs/` | Conda environments; each `envs/<name>.yaml` creates the environment `<name>` and has a `<name>.linux-64.pin.txt` lock |
| `reports/` | `report.html` (self-contained HTML report, from the input reads to the last step) and `sessionInfo.txt` |

Large files (FASTQ, primer-trimmed and filtered reads, SILVA) are never stored in the repository: they live under `$DATA`.

### Shared scripts

```
scripts/
├── html_report.py          # HTML report builder (same file as in bulk-rna-seq and proteomics)
└── amplicon/               # 16S amplicon projects
    ├── trim_primers.sh     # cutadapt, both mate orientations
    ├── quality_profiles.R, filter_trim.R, dada2_denoise.R, assign_taxonomy.R, build_phyloseq.R
    ├── analysis.R          # composition, alpha/beta diversity, PERMANOVA, DESeq2 + ANCOM-BC2
    │                       # (groups, comparisons and pairing come from each project's config.yaml)
    └── report.py           # HTML report of an amplicon project
```

## Requirements

- Linux (x86-64), Conda and Snakemake ≥ 9. Each rule creates its own environment with `--use-conda`; exact
  versions are in `envs/*.linux-64.pin.txt`.
- Internet access (ENA, Zenodo for SILVA).
- The `DATA` environment variable, pointing to storage for the FASTQ files, intermediate reads and SILVA.

```bash
git clone git@github.com:arlindomatias/microbiome.git
cd microbiome/ptx_rats_emtab13426      # or ptx_patients_emtab13423
export DATA=/path/to/large/storage
snakemake -n --use-conda                                   # dry run
snakemake --use-conda --cores 8 --resources downloads=4    # full run
```

## Data and citation

The data are public in ArrayExpress / ENA. Please cite the original authors when using them:
Zhong et al. (2023), *Blockade of CCR5 suppresses paclitaxel-induced peripheral neuropathic pain caused by
increased deoxycholic acid*, *Cell Reports* 42:113386.

## License

[MIT](LICENSE) © 2026 Arlindo César Matias Pereira. The license covers the code and the analysis outputs of this
repository; the input data (`data/`) belong to the original authors.
