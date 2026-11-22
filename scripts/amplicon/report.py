"""HTML report for a 16S amplicon project: reads → ASVs → composition → diversity → differential abundance.

Run through Snakemake (`script:` directive, rule `report`). Project-specific text comes from
`snakemake.params` (design, notes, pipeline, read_notes, error_notes, paper_taxa_text, label_notes) and the
study design (groups, comparisons, subject column) and cutoffs from the config. It only reads results/, plots/
and reports/sessionInfo.txt; every number in the text is computed from those tables.
"""
import os
import sys

import pandas as pd

# Shared report module: <repository>/scripts/html_report.py (this script lives in scripts/amplicon/).
sys.path.insert(0, os.path.join(snakemake.scriptdir, ".."))  # noqa: F821
from html_report import Report  # noqa: E402

inp, out, prm, cfg = snakemake.input, snakemake.output, snakemake.params, snakemake.config  # noqa: F821
P, R = "plots", "results"
q_cut = cfg["da_qvalue"]

GROUP_LABELS = dict(cfg["groups"])
REF_GROUP = next(iter(GROUP_LABELS))
COMPARISONS = [c["name"] for c in cfg["comparisons"]]
PAIRED = [c["name"] for c in cfg["comparisons"] if c.get("paired")]
UNPAIRED = [c for c in COMPARISONS if c not in PAIRED]
SUBJECT = cfg["subject_column"]
RUN_LABELS = dict(cfg.get("seq_run_labels", {}))


def param(name, default=""):
    """Optional project text from the rule's params."""
    return getattr(prm, name, default) or default


def names(items):
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


# How the two-group tests and models treat paired comparisons, for the text.
if PAIRED and UNPAIRED:
    WILCOXON_TEXT = (f"the rank-sum test for {names(UNPAIRED)}, and the signed-rank test, paired by {SUBJECT}, "
                     f"for {names(PAIRED)}")
elif PAIRED:
    WILCOXON_TEXT = f"the signed-rank test, paired by {SUBJECT}"
else:
    WILCOXON_TEXT = "the rank-sum test; all comparisons are between independent samples"


def read(name):
    path = os.path.join(R, f"{name}.tsv")
    try:
        return pd.read_csv(path, sep="\t")
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def fmt_p(p, prefix="p"):
    if pd.isna(p):
        return ""
    return f"{prefix} < 0.001" if p < 0.001 else f"{prefix} = {p:.2g}"


samples_all = read("samples")
EXCLUDED = dict(cfg.get("exclude_samples") or {})
samples = samples_all[~samples_all["sample"].isin(EXCLUDED)].reset_index(drop=True)
track = read("read_tracking")
summ = read("run_summary").set_index("key")["value"]
phylum_tests = read("phylum_tests")
alpha_tests = read("alpha_diversity_tests")
permanova = read("permanova")
d2c = read("distance_to_control")
d2c_tests = read("distance_to_control_tests")
comp = read("composition_group_means")
da = read("differential_abundance_genus")
paper = read("paper_taxa_tests")
label_check = read("authors_label_check")

n_ctrl = int((samples["group"] == REF_GROUP).sum())
ref_label = GROUP_LABELS[REF_GROUP].lower()

both = da[da["called_by"] == "both"]
robust_genera = both["genus"].nunique()
n_sig_comparisons = permanova[permanova["p_adj"] < 0.05].shape[0]

r = Report(prm.title, subtitle=prm.subtitle, code_files=list(inp.code))
r.kpis([
    (len(samples), "samples"),
    (f"{int(track['nonchim'].sum()):,}", "non-chimeric read pairs"),
    (f"{int(summ['asvs_kept']):,}", "bacterial ASVs"),
    (f"{int(summ['genera']):,}", "genera"),
    (f"{n_sig_comparisons} / {len(permanova)}", "PERMANOVA tests with padj < 0.05"),
    (f"{robust_genera}", "genera called by both DA methods"),
])

# ---------------------------------------------------------------- design
r.h2("Study design")
r.md(prm.design)
r.md(param("groups_text"))
cmp_lines = "\n".join(
    f"* **{c['name']}**: {GROUP_LABELS[c['groups'][1]]} vs {GROUP_LABELS[c['groups'][0]]}"
    + (f", paired by {SUBJECT} (grey lines in the plots join the same {SUBJECT})" if c.get("paired") else "")
    + "." for c in cfg["comparisons"])
r.md(f"""
Comparisons (fold changes are the second group relative to the first):

{cmp_lines}
""")
design = (samples.assign(group=pd.Categorical(samples["group"].map(GROUP_LABELS),
                                              categories=list(GROUP_LABELS.values())),
                         seq_run=samples["seq_run"].map(RUN_LABELS).fillna(samples["seq_run"]))
          .groupby(["group", "seq_run"], observed=True)
          .agg(samples=("sample", "size"), **{f"{SUBJECT}s": (SUBJECT, "nunique")}).reset_index()
          .rename(columns={"seq_run": "sequencing run"}))
r.table(design, "Samples per group and sequencing run.")
r.md(param("notes"))

r.h2("Pipeline")
r.md(prm.pipeline)
r.md("""
```
ENA FASTQ ──► cutadapt (both mate orientations) ──► DADA2 filterAndTrim
          ──► DADA2 learnErrors/dada per run × orientation batch ──► mergePairs ──► removeBimeraDenovo
          ──► SILVA 138.2 assignTaxonomy + addSpecies
          ──► MAFFT + FastTree (midpoint-rooted) ──► phyloseq ──► statistics ──► this report
```
""")

# ---------------------------------------------------------------- input and read processing
r.h2("Input data")
r.md(f"""
The input is {len(samples_all)} paired-end samples from ENA ({len(samples)} analysed). Sample names, groups, {SUBJECT}s and sequencing runs come
from the ArrayExpress sample sheet (`{cfg["sdrf"]}`); the full sheet is `results/samples.tsv`.
""")
if EXCLUDED:
    r.md("\n".join(f"> **{k} is excluded** from the analysis: {v}. It is kept in the read tracking table."
                    for k, v in EXCLUDED.items()))
r.table(samples_all, "Samples (including excluded ones).", collapsed=True)

r.h2("Read processing")
track_all = track
track = track[~track["sample"].isin(EXCLUDED)]  # read-processing numbers of the analysed samples
orient = track["orientation"].value_counts()
primer_pct = 100 * track["primers_found"] / track["raw"]
chim_pct = 100 * (1 - track["nonchim"].sum() / track["merged"].sum())
n_swapped = int(orient.get("swapped", 0))
orient_text = (f"in {n_swapped} of {len(track)} samples R1 starts with the reverse primer (806R)" if n_swapped
               else f"R1 starts with the forward primer (338F) in all {len(track)} samples")
r.md(f"""
* **Mate orientation:** {orient_text}. `trim_primers.sh` tries both orientations and re-orients the reads, so
  "fwd" reads always start after 338F. Both primers were found in {primer_pct.min():.1f}–{primer_pct.max():.1f}% of pairs per sample.
{param("read_notes").strip()}
* **Chimeras:** {chim_pct:.0f}% of merged reads were removed as chimeras. Between
  {int(track["nonchim"].min()):,} and {int(track["nonchim"].max()):,} non-chimeric pairs per sample remain
  ({track["retained_pct"].min():.0f}–{track["retained_pct"].max():.0f}% of the raw pairs).
""")
r.figure(f"{P}/read_tracking.png", "Read pairs per sample after each step.", width="85%")
r.table(track_all, "Read tracking (`results/read_tracking.tsv`).", collapsed=True)

r.h3("Quality profiles and error models")
r.md(param("error_notes",
           "Aggregate quality profiles of the primer-trimmed reads were used to choose the truncation lengths."))
r.figure(f"{P}/quality_profiles.png", "Aggregate DADA2 quality profiles of the primer-trimmed reads.", code=None)
r.figure(f"{P}/error_models.png", "DADA2 error models per sequencing run × orientation batch (points: observed; "
         "black line: fitted; red line: expected from the nominal quality).", code=None)

r.h3("ASV table and filtering")
r.md(f"""
DADA2 produced **{int(summ["asvs_total"]):,} ASVs**. {int(summ["asvs_nonbacterial"]):,} non-bacterial (Archaea or
unclassified kingdom) and {int(summ["asvs_organelle"]):,} chloroplast/mitochondria ASVs were removed, leaving
**{int(summ["asvs_kept"]):,}**. {summ["genus_assigned_pct"]}% of the kept ASVs ({summ["genus_assigned_reads_pct"]}%
of the reads) have a genus assignment. Libraries have {int(summ["reads_min"]):,}–{int(summ["reads_max"]):,} reads
(median {int(float(summ["reads_median"])):,}); for alpha diversity, Jaccard and unweighted UniFrac they are
rarefied to **{int(summ["rarefaction_depth"]):,} reads**.
""")
r.figure(f"{P}/rarefaction.png", "Rarefaction curves. Dashed line: rarefaction depth.", width="85%")

# ---------------------------------------------------------------- composition
r.h2("Taxonomic composition")
r.h3("Group averages")
r.figure(f"{P}/composition_group_means.png", "Mean relative abundance per group: phyla and the 8 most abundant "
         "genera (others pooled as Other).")
top_phyla = (comp[(comp["rank"] == "Phylum") & ~comp["taxon"].isin(["Other", "Unassigned"])]
             .pivot(index="taxon", columns="group", values="mean_abundance")
             .rename(columns=GROUP_LABELS).mul(100))
top_phyla = top_phyla.loc[top_phyla.mean(axis=1).sort_values(ascending=False).index]
top_phyla = top_phyla.apply(lambda c: c.map(lambda v: f"{v:.1f}")).reset_index()
r.table(top_phyla.rename(columns={"taxon": "phylum"}), "Mean relative abundance (%) of the main phyla.")

r.h3("Per sample")
r.figure(f"{P}/composition_per_sample.png", "Relative abundance per sample.")

r.h3("Main phyla and the Bacillota/Bacteroidota ratio")
r.md(f"""
Each panel shows one phylum (or the log2 Bacillota/Bacteroidota ratio). Brackets show BH-adjusted Wilcoxon
p-values: {WILCOXON_TEXT}.
""")
r.figure(f"{P}/phyla_boxplots.png", "Main phyla (relative abundance, %) and the log2 Bacillota/Bacteroidota ratio.")
r.table(phylum_tests.pivot(index="variable", columns="comparison", values="p_adj")[COMPARISONS].reset_index(),
        "BH-adjusted Wilcoxon p-values (`results/phylum_tests.tsv`).")

# ---------------------------------------------------------------- alpha
r.h2("Alpha diversity")
r.md(f"""
Samples are rarefied to {int(summ["rarefaction_depth"]):,} reads. Brackets show BH-adjusted Wilcoxon p-values
({WILCOXON_TEXT}).
""")
r.figure(f"{P}/alpha_diversity.png", "Observed ASVs, Chao1, Shannon and inverse Simpson indices.")
r.table(alpha_tests.pivot(index="index", columns="comparison", values="p_adj")[COMPARISONS].reset_index(),
        "BH-adjusted Wilcoxon p-values (`results/alpha_diversity_tests.tsv`).")

# ---------------------------------------------------------------- beta
r.h2("Beta diversity")
r.md(f"""
Four distances:

* **Bray-Curtis**: relative abundances.
* **Jaccard**: presence/absence, rarefied.
* **Unweighted UniFrac**: rarefied.
* **Weighted UniFrac**: relative abundances.

PERMANOVA (`adonis2`, {cfg["permutations"]:,} permutations) is run for each comparison separately.{
f" Permutations for {names(PAIRED)} are restricted within {SUBJECT}s." if PAIRED else ""} `betadisper` tests
whether groups differ in spread.
""")
r.h3("Ordination (PCoA)")
r.md(f"""
Each point is a sample.{f" Grey segments join samples of the same {SUBJECT}." if PAIRED else ""} Large outlined
points are group centroids, and dashed ellipses are 68% normal-theory regions.
""")
r.figure(f"{P}/pcoa.png", "PCoA for each distance, with PERMANOVA R² and BH-adjusted p-values.")

r.h3("Effect size across distances")
r.md("""
PERMANOVA R² is the fraction of between-sample variation explained by the grouping. Asterisks mark comparisons
where `betadisper` finds unequal spread (p < 0.05). In those cases part of the effect may reflect dispersion rather
than a shift in location.
""")
r.figure(f"{P}/permanova_r2.png", "PERMANOVA R² and BH-adjusted p-values.", width="85%")
r.table(permanova, f"PERMANOVA and betadisper (`results/permanova.tsv`; `p_adj`: BH over all {len(permanova)} "
        "tests).")

r.h3("Distance from the control community")
d2c_mean = d2c.groupby("group")["bray_curtis_to_controls"].mean().rename(index=GROUP_LABELS)
r.md(f"""
For every sample, this is the mean Bray-Curtis dissimilarity to the {n_ctrl} {ref_label} samples (excluding itself).
Values near 1 mean almost nothing is shared. Group means: {"; ".join(f"{g} {v:.2f}" for g, v in d2c_mean.items())}.
""")
r.figure(f"{P}/distance_to_control.png", f"Mean Bray-Curtis dissimilarity to the {ref_label} samples.",
         width="60%")

# ---------------------------------------------------------------- differential abundance
def model_text(unpaired, paired):
    if PAIRED and UNPAIRED:
        return f"`{unpaired}` for {names(UNPAIRED)}, `{paired}` for {names(PAIRED)}"
    return f"`{paired if PAIRED else unpaired}`"



r.h2("Differential abundance (genus level)")
r.md(f"""
ASVs are agglomerated to genus. Unassigned genera are labelled "&lt;lowest assigned rank&gt; (unclassified)". Each
comparison is fitted on its own subset, keeping genera present in at least {cfg["min_prevalence"]:.0%} of that
subset. Two methods are used:

* **DESeq2** (negative binomial, `poscounts`): {model_text("~ group", f"~ {SUBJECT} + group")}.
* **ANCOM-BC2** (bias-corrected log-linear model): {model_text("group", f"{SUBJECT} + group")}.

Genera with BH `q < {q_cut}` in **both** methods are treated as robust. Positive log2FC = higher in the second
group of the comparison.
""")
if isinstance(summ.get("ancombc2_failed"), str) and summ["ancombc2_failed"]:
    r.md(f"> ANCOM-BC2 failed for: {summ['ancombc2_failed']}. Those comparisons use DESeq2 only.")
counts = (da.groupby(["comparison", "called_by"]).size().unstack(fill_value=0)
          .reindex(index=COMPARISONS, columns=["both", "DESeq2 only", "ANCOM-BC2 only", "not significant"],
                   fill_value=0).reset_index())
counts.insert(1, "genera tested", counts.iloc[:, 1:].sum(axis=1))
r.table(counts, f"Genera per comparison by the method(s) that call them significant (q < {q_cut}).")

r.h3("Volcano plots")
r.md("""
Each point is a genus. The x axis is the DESeq2 fold change and the y axis its p-value. Colour shows which
method(s) call the genus significant. The six most significant genera called by both methods are labelled in each
panel.
""")
r.figure(f"{P}/volcano.png", "Volcano plots, one per comparison.")

r.h3("Genera that differ between treatments")
if robust_genera:
    r.md("""
These panels show the relative abundance (log scale) of every genus called significant by **both** methods in at
least one comparison. Brackets give the ANCOM-BC2 q-value for each comparison in which the genus was tested.
""")
    r.figure(f"{P}/robust_genera.png", "Genera significant in both DESeq2 and ANCOM-BC2.", optional=True)
    tbl = (both.sort_values(["comparison", "q_ancombc2"])
           [["comparison", "genus", "phylum", "log2FC_deseq2", "q_deseq2", "log2FC_ancombc2", "q_ancombc2",
             "passed_sensitivity"]]
           .assign(passed_sensitivity=lambda d: d["passed_sensitivity"].map({True: "yes", False: "no"})))
    r.table(tbl, "Genera called by both methods (`results/differential_abundance_genus.tsv` has every genus).")
else:
    r.md(f"No genus is called by both methods at `q < {q_cut}`.")

r.h3("Fold changes of significant genera")
r.figure(f"{P}/da_fold_changes.png", "log2 fold change of every genus significant in at least one method "
         "(ANCOM-BC2 when available, else DESeq2).", optional=True)

r.h3("Heatmap of the most abundant genera")
r.md("""
The 30 most abundant genera, each row scaled across samples (z-score of log10 relative abundance). Blue marks
below-average abundance for that genus and red above-average. Rows are clustered, and the bar on the left gives
the phylum.
""")
r.figure(f"{P}/heatmap_top_genera.png", "Top 30 genera, z-score of log10 relative abundance.")

# ---------------------------------------------------------------- paper taxa
r.h2("Taxa highlighted by the original study")
r.md(param("paper_taxa_text"))
r.md("""
All SILVA genera whose names contain one of these terms (and that are present in at least 25% of samples) are
shown, with BH-adjusted Wilcoxon p-values as above.
""")
r.figure(f"{P}/paper_taxa.png", "Genera highlighted by the original study (relative abundance, %).")
if len(paper):
    r.table(paper.pivot(index="genus", columns="comparison", values="p_adj")[COMPARISONS].reset_index(),
            "BH-adjusted Wilcoxon p-values (`results/paper_taxa_tests.tsv`).")

# ---------------------------------------------------------------- authors' table
r.h2("Consistency with the authors' processed table")
r.md(f"""
The authors' genus table (Greengenes, relative abundance; `{cfg["authors_genus_table"]}`) is compared with
this analysis for genera with the same name in both taxonomies. Unclassified groups and bracketed Greengenes names
(e.g. *[Ruminococcus]*, a Lachnospiraceae genus unrelated to SILVA *Ruminococcus*) are excluded, and SILVA suffixes
are stripped before matching. **{int(float(summ["authors_genera_matched"]))} genera** match, and their mean abundances correlate
with Spearman rho = **{float(summ["authors_rho_means"]):.2f}**.
""")
r.figure(f"{P}/authors_comparison.png", "Mean relative abundance per genus: authors' table vs this analysis.",
         width="70%")

r.h3("Sample label check")
same = label_check["best_authors_sample"] == label_check["sample"]
lc = label_check.merge(samples[["sample", "group", SUBJECT]], on="sample", how="left")
subject_of = dict(zip(samples["sample"], samples[SUBJECT]))
other_tp = (~same) & (lc["best_authors_sample"].map(subject_of) == lc[SUBJECT])
per_group = "; ".join(f"{GROUP_LABELS[g]} {int(same[lc['group'] == g].sum())} of {int((lc['group'] == g).sum())}"
                      for g in GROUP_LABELS)
other_txt = (f" In **{int(other_tp.sum())}** samples the best match is the **same {SUBJECT} under another label**."
             if other_tp.any() else "")
r.md(f"""
Each sample here is correlated (Spearman, shared genera) with every sample of the authors' table.

> **{int(same.sum())} of {len(label_check)} samples** best match the authors' sample with the same label
> ({per_group}).{other_txt} {param("label_notes").strip()}
""")
r.figure(f"{P}/label_check.png", "Spearman correlation between every sample of this analysis (rows) and of the "
         "authors' table (columns).")
r.table(label_check, "Best-matching authors' sample for each sample (`results/authors_label_check.tsv`).",
        collapsed=True)

# ---------------------------------------------------------------- appendix
r.appendix(
    outputs=[
        ("results/samples.tsv", "Sample sheet parsed from the SDRF"),
        ("results/read_tracking.tsv", "Read pairs per sample after each step"),
        ("results/taxonomy.tsv", "SILVA 138.2 taxonomy and sequence of every ASV"),
        ("results/asv_counts.tsv", "ASV count table with taxonomy"),
        ("data/asvs.fasta, data/asvs.nwk", "ASV sequences and the midpoint-rooted tree"),
        ("data/ps.rds", "phyloseq object (not versioned, regenerated by `build_phyloseq`)"),
        ("results/composition_group_means.tsv", "Mean relative abundance per group (phylum, top genera)"),
        ("results/genus_relative_abundance.tsv", "Relative abundance of every genus per sample"),
        ("results/phylum_tests.tsv, results/alpha_diversity*.tsv", "Phylum and alpha diversity values and tests"),
        ("results/permanova.tsv", "PERMANOVA and betadisper"),
        ("results/distance_to_control*.tsv", "Bray-Curtis distance to controls and tests"),
        ("results/deseq2_genus.tsv, results/ancombc2_genus.tsv", "Differential abundance per method"),
        ("results/differential_abundance_genus.tsv", "Both methods merged, with the call per genus"),
        ("results/paper_taxa_tests.tsv", "Tests for the genera highlighted by the original study"),
        ("results/authors_comparison.tsv, results/authors_label_check.tsv", "Comparison with the authors' table"),
        ("results/run_summary.tsv", "Key numbers of the run (read by this report)"),
    ],
    env_files=sorted(inp.envs),
    config_files=list(inp.config),
    session_file=inp.session,
)

r.write(out.html)
