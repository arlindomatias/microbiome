# Downstream statistics and figures for a 16S amplicon project (phyloseq object -> tables + plots).
# Composition, alpha diversity, beta diversity (PERMANOVA, betadisper), differential abundance
# at genus level (DESeq2 + ANCOM-BC2), taxa highlighted by the original study and a comparison
# with the authors' processed genus table.
#
# Run by the `analysis` rule (Snakemake `script:`). Writes result tables to results/, figures to
# plots/ and the R session to reports/sessionInfo.txt. The HTML report (report.py) only reads
# these files. Blocks between `# <figure name="...">` and `# </figure>` are shown as code in
# the report next to the figure of the same name.
#
# The study design comes from the config (through the rule's params): `groups` (ordered group
# keys and labels; the first group is the reference for "distance to control"), `comparisons`
# (name, two groups, paired or not), `subject` (metadata column identifying the individual;
# paired comparisons match samples on it) and `paper_taxa` (regex of genera highlighted by the
# original study).
log <-file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(permute)
  library(DESeq2)
  library(ANCOMBC)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

seed <- snakemake@params[["seed"]]
set.seed(seed)
threads <- snakemake@threads
results_dir <- snakemake@params[["results_dir"]]
plots_dir <- snakemake@params[["plots_dir"]]
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
da_q <- snakemake@params[["da_qvalue"]]
min_prev <- snakemake@params[["min_prevalence"]]
n_perm <- snakemake@params[["permutations"]]
subject_col <- snakemake@params[["subject"]]
paper_taxa <- snakemake@params[["paper_taxa"]]

# Tables keep the project's own name for the subject column (e.g. "rat", "patient").
save_table <- function(df, name) {
  names(df)[names(df) == "subject"] <- subject_col
  write_tsv(df, file.path(results_dir, paste0(name, ".tsv")))
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(plots_dir, paste0(name, ".png")), p, width = width, height = height,
         dpi = 120, bg = "white")
}
summary_rows <- list()
add_summary <- function(key, value) summary_rows[[key]] <<- as.character(value)

# Eight categorical slots of the validated palette, plus a neutral grey for "Other".
taxon_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
other_col <- "#b5b4ae"
# Group colours: the first slots of the same palette, in group order
# (all pairs pass colour-vision-deficiency checks).
group_labels <- unlist(snakemake@params[["groups"]])
group_levels <- names(group_labels)
ref_group <- group_levels[1]
group_cols <- setNames(taxon_cols[seq_along(group_levels)], group_levels)
scale_group_colour <- scale_colour_manual(values = group_cols, labels = group_labels, name = NULL)
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92"),
                  legend.position = "top"))

# ---------------------------------------------------------------- input and filtering
ps_raw <- readRDS(snakemake@input[["ps"]])
sample_data(ps_raw)$subject <- as.character(sample_data(ps_raw)[[subject_col]])
stopifnot(identical(levels(sample_data(ps_raw)$group), group_levels))
meta <-data.frame(sample_data(ps_raw)) %>% rownames_to_column("sample_id")

# Keep bacteria only (drop Archaea, unclassified kingdom, chloroplast, mitochondria).
tax_raw <- as.data.frame(tax_table(ps_raw))
is_bact <- !is.na(tax_raw$Kingdom) & tax_raw$Kingdom == "Bacteria"
is_organelle <- (!is.na(tax_raw$Order) & tax_raw$Order == "Chloroplast") |
  (!is.na(tax_raw$Family) & tax_raw$Family == "Mitochondria")
ps <- prune_taxa(taxa_names(ps_raw)[is_bact & !is_organelle], ps_raw)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
depth <- min(sample_sums(ps))
ps_rare <- rarefy_even_depth(ps, sample.size = depth, rngseed = seed,
                             replace = FALSE, verbose = FALSE)

add_summary("asvs_total", ntaxa(ps_raw))
add_summary("asvs_nonbacterial", sum(!is_bact))
add_summary("asvs_organelle", sum(is_bact & is_organelle))
add_summary("asvs_kept", ntaxa(ps))
add_summary("reads_min", min(sample_sums(ps)))
add_summary("reads_max", max(sample_sums(ps)))
add_summary("reads_median", median(sample_sums(ps)))
add_summary("rarefaction_depth", depth)
add_summary("genus_assigned_pct", round(100 * mean(!is.na(tax_table(ps)[, "Genus"])), 1))
add_summary("genus_assigned_reads_pct",
            round(100 * sum(taxa_sums(ps)[!is.na(tax_table(ps)[, "Genus"])]) / sum(taxa_sums(ps)), 1))

# Comparisons from the config: list(name, groups = c(first, second), paired). Fold changes are
# second vs first group.
comparisons <- lapply(snakemake@params[["comparisons"]], function(x)
  list(groups = unlist(x[["groups"]]), paired = isTRUE(x[["paired"]])))
names(comparisons) <- vapply(snakemake@params[["comparisons"]], function(x) x[["name"]], "")
stopifnot(all(unlist(lapply(comparisons, `[[`, "groups")) %in% group_levels))
# Subjects sampled more than once (joined by grey lines in the plots).
has_pairs <- any(duplicated(meta$subject))
# Samples ordered by group, then by the number in the subject name (C2 before C10).
sample_order <- meta %>% arrange(group, as.integer(gsub("\\D", "", subject)), sample_id) %>% pull(sample_id)
n_samples <- nrow(meta)

# One Wilcoxon test per comparison and variable: rank-sum for unpaired comparisons,
# signed-rank (paired by subject) for paired ones; BH-adjusted within the variable.
group_tests <- function(d) {
  p <- vapply(comparisons, function(cmp) {
    a <- d %>% filter(group == cmp$groups[1]) %>% arrange(subject)
    b <- d %>% filter(group == cmp$groups[2]) %>% arrange(subject)
    if (cmp$paired) {
      stopifnot(identical(a$subject, b$subject))
      wilcox.test(a$value, b$value, paired = TRUE, exact = FALSE)$p.value
    } else {
      wilcox.test(a$value, b$value, exact = FALSE)$p.value
    }
  }, numeric(1))
  tibble(comparison = names(comparisons), p = unname(p)) %>% mutate(p_adj = p.adjust(p, "BH"))
}

fmt_p <- function(p, prefix = "p") {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, paste0(prefix, " < 0.001"),
                paste0(prefix, " = ", trimws(formatC(signif(p, 2), format = "fg", digits = 2)))))
}

# Significance brackets on a group box plot (x = groups in config order); short brackets are
# drawn first (lowest), then longer ones above them.
# stats: columns facet, comparison, label; values: columns facet, value (plotted scale).
bracket_pos <- tibble(
  comparison = names(comparisons),
  a = vapply(comparisons, function(cmp) match(cmp$groups[1], group_levels), 1L),
  b = vapply(comparisons, function(cmp) match(cmp$groups[2], group_levels), 1L)
) %>%
  transmute(comparison, x = pmin(a, b), xend = pmax(a, b)) %>%
  arrange(xend - x, x) %>% mutate(level = row_number())
# Grey lines joining samples of the same subject (paired designs only).
subject_lines <- function(geom, colour) {
  if (!has_pairs) return(NULL)
  geom(data = ~ .x %>% group_by(across(any_of("facet")), subject) %>% filter(n() > 1) %>%
         ungroup() %>% arrange(subject, group),
       aes(group = subject), colour = colour, linewidth = 0.4)
}
# Without repeated subjects, points are jittered so that large groups stay readable.
point_position <- if (has_pairs) "identity" else position_jitter(width = 0.15, height = 0, seed = 1)
add_brackets <- function(stats, values) {
  rng <- values %>% group_by(facet) %>%
    summarise(top = max(value, na.rm = TRUE),
              span = max(diff(range(value, na.rm = TRUE)), 1e-9), .groups = "drop")
  b <- stats %>% inner_join(bracket_pos, by = "comparison") %>%
    left_join(rng, by = "facet") %>%
    mutate(y = top + span * (0.07 + 0.12 * (level - 1)))
  list(
    geom_segment(data = b, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE,
                 colour = "grey45", linewidth = 0.35),
    geom_text(data = b, aes(x = (x + xend) / 2, y = y, label = label), inherit.aes = FALSE,
              vjust = -0.35, size = 3, colour = "grey20")
  )
}

# Box plot of a per-sample value by group, samples of the same subject joined.
group_boxplot <- function(d, stats, ncol = 4, ylab = NULL) {
  ggplot(d, aes(group, value, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.5) +
    subject_lines(geom_line, "grey78") +
    geom_point(size = 2, position = point_position, alpha = if (has_pairs) 1 else 0.75) +
    add_brackets(stats, d) +
    facet_wrap(~ facet, scales = "free_y", ncol = ncol, labeller = label_wrap_gen(width = 26)) +
    scale_group_colour +
    scale_x_discrete(labels = group_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = NULL, y = ylab) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
}

# ---------------------------------------------------------------- read tracking and depth
# Samples excluded in build_phyloseq (config exclude_samples) are left out of the plot.
track <- read_tsv(snakemake@input[["track"]], show_col_types = FALSE) %>%
  filter(sample %in% meta$sample_id)

# <figure name="read_tracking">
p_track <- track %>% left_join(meta %>% select(sample = sample_id, group), by = "sample") %>%
  select(sample, group, raw, primers_found, filtered, merged, nonchim) %>%
  pivot_longer(raw:nonchim, names_to = "step", values_to = "reads") %>%
  mutate(step = factor(step, levels = c("raw", "primers_found", "filtered", "merged", "nonchim"))) %>%
  ggplot(aes(step, reads, group = sample, colour = group)) +
  geom_line(alpha = 0.6) + geom_point(size = 1.4) +
  scale_group_colour + scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  labs(x = NULL, y = "Read pairs")
save_plot(p_track, "read_tracking", 9, 4.5)
# </figure>

# <figure name="rarefaction">
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
p_rare <- rarecurve(otu, step = 500, tidy = TRUE) %>%
  left_join(meta, by = c(Site = "sample_id")) %>%
  ggplot(aes(Sample, Species, group = Site, colour = group)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(xintercept = depth, linetype = "dashed", colour = "grey50") +
  scale_group_colour +
  labs(x = "Reads sampled", y = "Observed ASVs", caption = "Dashed line: rarefaction depth")
save_plot(p_rare, "rarefaction", 9, 4.5)
# </figure>

# ---------------------------------------------------------------- composition
# Relative abundance per sample at a rank: top-n taxa by mean abundance, the rest pooled as
# "Other"; taxa without an assignment at that rank are labelled "Unassigned".
composition <- function(ps, rank, top_n = 8) {
  df <- psmelt(tax_glom(transform_sample_counts(ps, function(x) x / sum(x)), rank, NArm = FALSE))
  df[[rank]] <- ifelse(is.na(df[[rank]]), "Unassigned", df[[rank]])
  df <- df %>% group_by(Sample, group, subject, taxon = .data[[rank]]) %>%
    summarise(abundance = sum(Abundance), .groups = "drop")
  top <- df %>% filter(taxon != "Unassigned") %>% group_by(taxon) %>%
    summarise(m = mean(abundance)) %>% slice_max(m, n = top_n) %>% pull(taxon)
  df %>% mutate(taxon = ifelse(taxon %in% top, taxon, "Other")) %>%
    group_by(Sample, group, subject, taxon) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    mutate(taxon = factor(taxon, levels = c(top, "Other")))
}
taxon_fill <- function(lv) {
  scale_fill_manual(values = setNames(c(taxon_cols[seq_len(length(lv) - 1)], other_col), lv),
                    name = NULL)
}

phylum_df <- composition(ps, "Phylum")
genus_df <- composition(ps, "Genus")
save_table(bind_rows(phylum_df %>% mutate(rank = "Phylum"), genus_df %>% mutate(rank = "Genus")) %>%
             group_by(rank, group, taxon) %>% summarise(mean_abundance = mean(abundance), .groups = "drop") %>%
             arrange(rank, group, desc(mean_abundance)),
           "composition_group_means")

# <figure name="composition_group_means">
plot_group_means <- function(df, title) {
  m <- df %>% group_by(group, taxon) %>% summarise(abundance = mean(abundance), .groups = "drop")
  # Stack positions are computed over all segments before hiding labels of small ones.
  labs_df <- m %>% arrange(group, desc(taxon)) %>%
    group_by(group) %>% mutate(y = cumsum(abundance) - abundance / 2) %>% ungroup() %>%
    filter(abundance >= 0.04)
  ggplot(m, aes(group, abundance, fill = taxon)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.4) +
    geom_text(data = labs_df, aes(y = y, label = scales::percent(abundance, accuracy = 1)),
              colour = "white", size = 3.2, fontface = "bold") +
    taxon_fill(levels(m$taxon)) +
    scale_x_discrete(labels = group_labels) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = "Mean relative abundance", title = title) +
    theme(panel.grid.major.x = element_blank(), legend.position = "right")
}
save_plot(plot_group_means(phylum_df, "Phylum") + plot_group_means(genus_df, "Genus (top 8)"),
          "composition_group_means", 12, 5.5)
# </figure>

# <figure name="composition_per_sample">
plot_per_sample <- function(df, title) {
  ggplot(df %>% mutate(Sample = factor(Sample, levels = sample_order)), aes(Sample, abundance, fill = taxon)) +
    geom_col(width = 0.9, colour = "white", linewidth = 0.2) +
    facet_grid(~ group, scales = "free_x", space = "free_x",
               labeller = as_labeller(group_labels)) +
    taxon_fill(levels(df$taxon)) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = "Relative abundance", title = title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.grid.major.x = element_blank(), legend.position = "right")
}
save_plot(plot_per_sample(phylum_df, "Phylum") / plot_per_sample(genus_df, "Genus (top 8)"),
          "composition_per_sample", max(10, 3 + 0.12 * n_samples), 9)
# </figure>

# <figure name="phyla_boxplots">
phyla <- phylum_df %>% filter(!taxon %in% c("Other", "Unassigned")) %>%
  group_by(taxon) %>% filter(mean(abundance) >= 0.005) %>% ungroup() %>%
  transmute(sample = Sample, group, subject, facet = as.character(taxon), value = 100 * abundance)
fb <- phylum_df %>% filter(taxon %in% c("Bacillota", "Bacteroidota")) %>%
  select(Sample, group, subject, taxon, abundance) %>%
  pivot_wider(names_from = taxon, values_from = abundance) %>%
  transmute(sample = Sample, group, subject, facet = "log2(Bacillota / Bacteroidota)",
            value = log2(Bacillota / Bacteroidota))
phyla_d <- bind_rows(phyla, fb) %>%
  mutate(facet = factor(facet, levels = c(sort(unique(phyla$facet)), "log2(Bacillota / Bacteroidota)")))
phyla_stats <- phyla_d %>% group_by(facet) %>% group_modify(~ group_tests(.x)) %>% ungroup() %>%
  mutate(label = fmt_p(p_adj))
save_plot(group_boxplot(phyla_d, phyla_stats, ncol = 3, ylab = "Relative abundance (%) / log2 ratio"),
          "phyla_boxplots", 10, 7.5)
# </figure>
save_table(phyla_stats %>% rename(variable = facet) %>% select(-label), "phylum_tests")

# ---------------------------------------------------------------- alpha diversity
# <figure name="alpha_diversity">
alpha <- estimate_richness(ps_rare, measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) %>%
  mutate(sample_id = sample_names(ps_rare)) %>%
  left_join(meta, by = "sample_id") %>%
  select(sample = sample_id, group, subject, Observed, Chao1, Shannon, InvSimpson) %>%
  pivot_longer(c(Observed, Chao1, Shannon, InvSimpson), names_to = "facet", values_to = "value") %>%
  mutate(facet = factor(facet, levels = c("Observed", "Chao1", "Shannon", "InvSimpson")))
alpha_stats <- alpha %>% group_by(facet) %>% group_modify(~ group_tests(.x)) %>% ungroup() %>%
  mutate(label = fmt_p(p_adj))
save_plot(group_boxplot(alpha, alpha_stats, ncol = 4), "alpha_diversity", 10, 5)
# </figure>
save_table(alpha %>% rename(index = facet), "alpha_diversity")
save_table(alpha_stats %>% rename(index = facet) %>% select(-label), "alpha_diversity_tests")

# ---------------------------------------------------------------- beta diversity
dists <- list(
  "Bray-Curtis" = phyloseq::distance(ps_rel, "bray"),
  "Jaccard" = phyloseq::distance(ps_rare, "jaccard", binary = TRUE),
  "Unweighted UniFrac" = phyloseq::UniFrac(ps_rare, weighted = FALSE),
  "Weighted UniFrac" = phyloseq::UniFrac(ps_rel, weighted = TRUE)
)

# PERMANOVA per comparison; permutations restricted within rats for the paired comparison.
permanova_one <- function(d, cmp) {
  m <- meta %>% filter(group %in% cmp$groups) %>% mutate(group = droplevels(group))
  dm <- as.dist(as.matrix(d)[m$sample_id, m$sample_id])
  set.seed(seed)
  fit <- if (cmp$paired) {
    adonis2(dm ~ group, data = m, permutations = how(nperm = n_perm, blocks = factor(m$subject)))
  } else {
    adonis2(dm ~ group, data = m, permutations = n_perm)
  }
  disp <- permutest(betadisper(dm, m$group), permutations = n_perm)
  tibble(R2 = fit$R2[1], F = fit$F[1], p = fit$`Pr(>F)`[1], dispersion_p = disp$tab$`Pr(>F)`[1])
}
permanova <- expand_grid(distance = names(dists), comparison = names(comparisons)) %>%
  rowwise() %>%
  mutate(res = list(permanova_one(dists[[distance]], comparisons[[comparison]]))) %>%
  unnest(res) %>% ungroup() %>%
  mutate(p_adj = p.adjust(p, "BH"))
save_table(permanova, "permanova")

# <figure name="pcoa">
pcoa_plot <- function(d, name) {
  ord <- cmdscale(d, k = 2, eig = TRUE)
  var_expl <- 100 * ord$eig[1:2] / sum(ord$eig[ord$eig > 0])
  df <- data.frame(ord$points, sample_id = rownames(ord$points)) %>% left_join(meta, by = "sample_id")
  cent <- df %>% group_by(group) %>% summarise(X1 = mean(X1), X2 = mean(X2))
  st <- permanova %>% filter(distance == name) %>%
    mutate(txt = sprintf("%s: R² = %.2f, %s", sub(" \\(paired\\)", "", comparison), R2, fmt_p(p_adj)))
  ggplot(df, aes(X1, X2, colour = group)) +
    subject_lines(geom_path, "grey82") +
    stat_ellipse(level = 0.68, linewidth = 0.5, linetype = "dashed") +
    geom_point(size = 2.4, alpha = 0.9) +
    geom_point(data = cent, size = 6, shape = 21, fill = group_cols[as.character(cent$group)],
               colour = "white", stroke = 1.5) +
    scale_group_colour +
    labs(title = name, subtitle = paste(st$txt, collapse = "\n"),
         x = sprintf("PCoA1 (%.1f%%)", var_expl[1]), y = sprintf("PCoA2 (%.1f%%)", var_expl[2])) +
    theme(plot.subtitle = element_text(size = 9, colour = "grey30"))
}
save_plot(wrap_plots(Map(pcoa_plot, dists, names(dists)), ncol = 2, guides = "collect") &
            theme(legend.position = "top"),
          "pcoa", 10, 10)
# </figure>

# <figure name="permanova_r2">
p_r2 <- permanova %>%
  mutate(label = sprintf("R² = %.2f\n%s%s", R2, fmt_p(p_adj), ifelse(dispersion_p < 0.05, " *", "")),
         comparison = factor(comparison, levels = names(comparisons)),
         distance = factor(distance, levels = rev(names(dists)))) %>%
  ggplot(aes(comparison, distance, fill = R2)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label, colour = R2 > 0.35), size = 3.3, lineheight = 0.9) +
  scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = c(0, NA), name = "R²") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), legend.position = "right")
save_plot(p_r2, "permanova_r2", 10, 3.8)
# </figure>

# <figure name="distance_to_control">
bc <- as.matrix(dists[["Bray-Curtis"]])
ctrl_ids <- meta$sample_id[meta$group == ref_group]
d2c <- meta %>% rowwise() %>%
  mutate(value = mean(bc[sample_id, setdiff(ctrl_ids, sample_id)])) %>% ungroup() %>%
  transmute(sample = sample_id, group, subject, facet = "Bray-Curtis to controls", value)
d2c_stats <- d2c %>% group_by(facet) %>% group_modify(~ group_tests(.x)) %>% ungroup() %>%
  mutate(label = fmt_p(p_adj))
save_plot(group_boxplot(d2c, d2c_stats, ncol = 1, ylab = "Mean Bray-Curtis dissimilarity"),
          "distance_to_control", 6, 4.5)
# </figure>
save_table(d2c %>% select(sample, group, subject, bray_curtis_to_controls = value), "distance_to_control")
save_table(d2c_stats %>% select(-facet, -label), "distance_to_control_tests")

# ---------------------------------------------------------------- differential abundance (genus)
ps_genus <- tax_glom(ps, "Genus", NArm = FALSE)
gtax <- as.data.frame(tax_table(ps_genus))
genus_label <- apply(gtax[, c("Phylum", "Class", "Order", "Family", "Genus")], 1, function(r) {
  if (!is.na(r["Genus"])) return(r["Genus"])
  known <- r[!is.na(r)]
  if (length(known) == 0) "Bacteria (unclassified)" else paste0(tail(known, 1), " (unclassified)")
})
taxa_names(ps_genus) <- make.unique(unname(genus_label))
genus_phylum <- setNames(gtax$Phylum, taxa_names(ps_genus))
add_summary("genera", ntaxa(ps_genus))

prevalence_filter <- function(ps) {
  counts <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) counts <- t(counts)
  prune_taxa(rowMeans(counts > 0) >= min_prev, ps)
}
subset_cmp <- function(cmp) {
  p <- prune_samples(meta$sample_id[meta$group %in% cmp$groups], ps_genus)
  sample_data(p)$group <- droplevels(sample_data(p)$group)
  sample_data(p)$subject <- factor(sample_data(p)$subject)
  prevalence_filter(p)
}

genus_rel <- psmelt(transform_sample_counts(ps_genus, function(x) x / sum(x))) %>%
  select(sample_id = Sample, genus = OTU, abundance = Abundance) %>%
  left_join(meta, by = "sample_id")
save_table(genus_rel %>% select(sample = sample_id, group, genus, relative_abundance = abundance) %>%
             mutate(phylum = genus_phylum[genus]),
           "genus_relative_abundance")

run_deseq2 <- function(cmp, name) {
  p <- subset_cmp(cmp)
  dds <- phyloseq_to_deseq2(p, if (cmp$paired) ~ subject + group else ~ group)
  dds <- DESeq(dds, sfType = "poscounts", fitType = "local", quiet = TRUE)
  res <- results(dds, name = resultsNames(dds)[length(resultsNames(dds))], alpha = da_q)
  as.data.frame(res) %>% rownames_to_column("genus") %>%
    transmute(comparison = name, genus, baseMean, log2FC = log2FoldChange, lfcSE,
              p = pvalue, q = padj)
}
deseq <- bind_rows(Map(run_deseq2, comparisons, names(comparisons)))
save_table(deseq, "deseq2_genus")

ancom_failed <- character(0)
run_ancombc2 <- function(cmp, name) {
  p <- subset_cmp(cmp)
  counts <- as(otu_table(p), "matrix")
  if (!taxa_are_rows(p)) counts <- t(counts)
  set.seed(seed)
  # Count matrix + metadata input (a phyloseq input would need the 'microbiome' package).
  # Paired design: subject as a fixed effect (a random intercept cannot be fitted with two
  # observations per subject).
  fit <- tryCatch(
    ancombc2(data = counts, taxa_are_rows = TRUE, meta_data = data.frame(sample_data(p)),
             fix_formula = if (cmp$paired) "subject + group" else "group",
             p_adj_method = "BH", prv_cut = 0, lib_cut = 0, group = "group",
             struc_zero = FALSE, neg_lb = FALSE, alpha = da_q, global = FALSE,
             n_cl = threads, verbose = FALSE),
    error = function(e) { cat(name, ": ANCOM-BC2 failed: ", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) {
    ancom_failed <<- c(ancom_failed, name)
    return(NULL)
  }
  r <- fit$res
  suffix <- sub("^lfc_", "", grep("^lfc_group", names(r), value = TRUE)[1])
  tibble(comparison = name, genus = r$taxon,
         # ANCOM-BC2 reports natural-log fold changes; converted to log2.
         log2FC = r[[paste0("lfc_", suffix)]] / log(2),
         se = r[[paste0("se_", suffix)]] / log(2),
         p = r[[paste0("p_", suffix)]], q = r[[paste0("q_", suffix)]],
         passed_sensitivity = r[[paste0("passed_ss_", suffix)]])
}
ancom <- bind_rows(Map(run_ancombc2, comparisons, names(comparisons)))
save_table(ancom, "ancombc2_genus")
add_summary("ancombc2_failed", paste(ancom_failed, collapse = "; "))

da <- full_join(
  deseq %>% select(comparison, genus, log2FC_deseq2 = log2FC, p_deseq2 = p, q_deseq2 = q),
  ancom %>% select(comparison, genus, log2FC_ancombc2 = log2FC, q_ancombc2 = q, passed_sensitivity),
  by = c("comparison", "genus")
) %>%
  mutate(phylum = genus_phylum[genus],
         sig_deseq2 = !is.na(q_deseq2) & q_deseq2 < da_q,
         sig_ancombc2 = !is.na(q_ancombc2) & q_ancombc2 < da_q,
         called_by = case_when(sig_deseq2 & sig_ancombc2 ~ "both",
                               sig_deseq2 ~ "DESeq2 only",
                               sig_ancombc2 ~ "ANCOM-BC2 only",
                               TRUE ~ "not significant"),
         called_by = factor(called_by, levels = c("both", "DESeq2 only", "ANCOM-BC2 only",
                                                  "not significant")),
         comparison = factor(comparison, levels = names(comparisons)))
save_table(da, "differential_abundance_genus")
called_cols <- c("both" = "#2a78d6", "DESeq2 only" = "#eb6834", "ANCOM-BC2 only" = "#1baf7a",
                 "not significant" = "#c9c8c2")

# <figure name="volcano">
p_volcano <- ggplot(da %>% filter(!is.na(p_deseq2)), aes(log2FC_deseq2, -log10(p_deseq2), colour = called_by)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(data = ~ filter(.x, called_by == "not significant"), size = 1.6, alpha = 0.7) +
  geom_point(data = ~ filter(.x, called_by != "not significant"), size = 2.4) +
  geom_text(data = ~ filter(.x, called_by == "both") %>% group_by(comparison) %>%
              slice_min(p_deseq2, n = 6) %>% ungroup(),
            aes(label = genus), size = 2.9, vjust = -0.8, colour = "grey15", check_overlap = TRUE) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_colour_manual(values = called_cols, name = sprintf("Significant (q < %s) in", da_q)) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  labs(x = "log2 fold change (DESeq2); positive = higher in the second group",
       y = "-log10 p (DESeq2)")
save_plot(p_volcano, "volcano", 13, 6)
# </figure>

robust <- da %>% filter(called_by == "both") %>% distinct(genus) %>% pull(genus)
if (length(robust) > 0) {
  # <figure name="robust_genera">
  pseudo <- 1e-3  # percent; zeros are plotted at this floor
  rob_d <- genus_rel %>% filter(genus %in% robust) %>%
    transmute(sample = sample_id, group, subject, facet = genus, value = log10(100 * abundance + pseudo))
  rob_stats <- da %>% filter(genus %in% robust, !is.na(q_ancombc2)) %>%
    transmute(facet = genus, comparison = as.character(comparison),
              label = paste0(fmt_p(q_ancombc2, "q"), ifelse(called_by == "both", " *", "")))
  # Order panels: largest ANCOM-BC2 fold change first.
  ord <- da %>% filter(genus %in% robust) %>% group_by(genus) %>%
    summarise(k = max(abs(log2FC_ancombc2), na.rm = TRUE)) %>% arrange(desc(k)) %>% pull(genus)
  rob_d$facet <- factor(rob_d$facet, levels = ord)
  rob_stats$facet <- factor(rob_stats$facet, levels = ord)
  n_rows <- ceiling(length(robust) / 5)
  p_rob <- group_boxplot(rob_d, rob_stats, ncol = 5, ylab = "Relative abundance (%)") +
    scale_y_continuous(breaks = log10(c(0, 0.01, 0.1, 1, 10, 100) + pseudo),
                       labels = c("0", "0.01", "0.1", "1", "10", "100"),
                       expand = expansion(mult = c(0.05, 0.22))) +
    labs(caption = "* significant in both DESeq2 and ANCOM-BC2. Zeros plotted at 0.001%.")
  save_plot(p_rob, "robust_genera", 12, 1 + 2.6 * n_rows)
  # </figure>
}

da_plot <- da %>% filter(called_by != "not significant") %>%
  mutate(log2FC = coalesce(log2FC_ancombc2, log2FC_deseq2))
if (nrow(da_plot) > 0) {
  # <figure name="da_fold_changes">
  p_da <- ggplot(da_plot, aes(log2FC, reorder(genus, log2FC), colour = called_by)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_segment(aes(x = 0, xend = log2FC, yend = reorder(genus, log2FC)), linewidth = 0.5) +
    geom_point(size = 2.4) +
    facet_grid(comparison ~ ., scales = "free_y", space = "free_y",
               labeller = label_wrap_gen(width = 14)) +
    scale_colour_manual(values = called_cols, name = "Significant in") +
    labs(x = "log2 fold change (ANCOM-BC2 when available, else DESeq2)\npositive = higher in the second group of the comparison",
         y = NULL) +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
  save_plot(p_da, "da_fold_changes", 10, max(4, 1.5 + 0.16 * nrow(da_plot)))
  # </figure>
}

# <figure name="heatmap_top_genera">
top30 <- genus_rel %>% group_by(genus) %>% summarise(m = mean(abundance)) %>%
  slice_max(m, n = 30) %>% pull(genus)
hm <- genus_rel %>% filter(genus %in% top30) %>%
  mutate(v = log10(abundance + 1e-5)) %>%
  group_by(genus) %>% mutate(z = (v - mean(v)) / sd(v)) %>% ungroup()
zmat <- hm %>% select(genus, sample_id, z) %>% pivot_wider(names_from = sample_id, values_from = z) %>%
  column_to_rownames("genus") %>% as.matrix()
row_order <- rownames(zmat)[hclust(dist(zmat), method = "average")$order]
hm <- hm %>% mutate(genus = factor(genus, levels = row_order),
                    sample_id = factor(sample_id, levels = sample_order),
                    z = pmax(pmin(z, 2.5), -2.5))
p_heat <- ggplot(hm, aes(sample_id, genus, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(~ group, scales = "free_x", space = "free_x", labeller = as_labeller(group_labels)) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
                       limits = c(-2.5, 2.5), name = "z-score") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_blank(), panel.grid = element_blank(), legend.position = "right")
phy_df <- tibble(genus = factor(row_order, levels = row_order), phylum = genus_phylum[row_order])
top_phy <- names(sort(table(phy_df$phylum), decreasing = TRUE))
phy_df$phylum <- factor(ifelse(phy_df$phylum %in% top_phy[1:min(4, length(top_phy))],
                               phy_df$phylum, "Other"),
                        levels = c(top_phy[1:min(4, length(top_phy))], "Other"))
p_phy <- ggplot(phy_df, aes(1, genus, fill = phylum)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = setNames(c(taxon_cols[seq_len(nlevels(phy_df$phylum) - 1)], other_col),
                                      levels(phy_df$phylum)), name = "Phylum") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), legend.position = "left")
save_plot(p_phy + p_heat + plot_layout(widths = c(0.04, 1), guides = "collect"),
          "heatmap_top_genera", max(12, 4 + 0.12 * n_samples), 9)
# </figure>

# ---------------------------------------------------------------- taxa highlighted by the study
# <figure name="paper_taxa">
focus <- genus_rel %>% filter(grepl(paper_taxa, genus)) %>%
  group_by(genus) %>% filter(mean(abundance > 0) >= 0.25) %>% ungroup() %>%
  transmute(sample = sample_id, group, subject, facet = genus, value = 100 * abundance)
focus_stats <- focus %>% group_by(facet) %>% group_modify(~ group_tests(.x)) %>% ungroup() %>%
  mutate(label = fmt_p(p_adj))
save_plot(group_boxplot(focus, focus_stats, ncol = 3, ylab = "Relative abundance (%)"),
          "paper_taxa", 10, 1.5 + 3 * ceiling(n_distinct(focus$facet) / 3))
# </figure>
save_table(focus_stats %>% rename(genus = facet) %>% select(-label), "paper_taxa_tests")

# ---------------------------------------------------------------- authors' processed table
# First column: taxon (full Greengenes lineage or genus name); one column per sample with
# relative abundances. The file may be UTF-16 (byte-order mark), as exported by Excel.
read_authors <- function(path) {
  bom <- readBin(path, "raw", 2)
  enc <- if (identical(bom, as.raw(c(0xff, 0xfe))) || identical(bom, as.raw(c(0xfe, 0xff)))) "UTF-16" else ""
  d <- read.delim(path, check.names = FALSE, fileEncoding = enc, blank.lines.skip = TRUE)
  names(d)[1] <- "Taxon"
  d %>% filter(!is.na(Taxon), Taxon != "")
}
# Greengenes names in brackets ([Ruminococcus], [Prevotella], [Eubacterium]) are taxa placed in
# another family by later taxonomies, so they are not matched to the SILVA genus of the same
# name. Greengenes "Shigella" corresponds to SILVA "Escherichia-Shigella".
authors_synonyms <- c(Shigella = "Escherichia-Shigella", Escherichia = "Escherichia-Shigella")
authors <- read_authors(snakemake@input[["authors_genus"]])
authors_long <- authors %>%
  mutate(genus = sub(".*g__", "", Taxon)) %>%
  filter(!grepl("^Unclassified|^$|^\\[", genus)) %>%
  mutate(genus = coalesce(authors_synonyms[genus], genus)) %>%
  select(-Taxon) %>%
  pivot_longer(-genus, names_to = "sample_id", values_to = "authors") %>%
  mutate(sample_id = gsub(".", "_", sample_id, fixed = TRUE)) %>%
  group_by(sample_id, genus) %>% summarise(authors = sum(authors), .groups = "drop")
ours_long <- genus_rel %>%
  filter(!grepl("unclassified", genus)) %>%
  mutate(genus = sub("\\.\\d+$", "", genus),
         genus = sub("_sensu_stricto.*|_group$|_\\d+$", "", genus),
         genus = gsub("\\[|\\]", "", genus)) %>%
  group_by(sample_id, genus) %>% summarise(ours = sum(abundance), .groups = "drop")
both <- inner_join(authors_long, ours_long, by = c("sample_id", "genus"))
mean_ab <- both %>% group_by(genus) %>% summarise(authors = mean(authors), ours = mean(ours))
rho_means <- cor(mean_ab$authors, mean_ab$ours, method = "spearman")
add_summary("authors_genera_matched", nrow(mean_ab))
add_summary("authors_rho_means", round(rho_means, 2))
save_table(mean_ab %>% rename(mean_authors = authors, mean_this_analysis = ours), "authors_comparison")

# <figure name="authors_comparison">
p_auth <- ggplot(mean_ab, aes(authors + 1e-5, ours + 1e-5)) +
  geom_abline(colour = "grey60", linetype = "dashed") +
  geom_point(colour = "#2a78d6", size = 2) +
  geom_text(data = slice_max(mean_ab, pmax(authors, ours), n = 10), aes(label = genus),
            size = 3, vjust = -0.8, colour = "grey25", check_overlap = TRUE) +
  scale_x_log10(labels = scales::percent) + scale_y_log10(labels = scales::percent) +
  labs(x = "Mean relative abundance, authors (Greengenes)",
       y = "Mean relative abundance, this analysis (SILVA)",
       caption = sprintf("%d genera matched; Spearman rho of means = %.2f; dashed line: identity",
                         nrow(mean_ab), rho_means))
save_plot(p_auth, "authors_comparison", 7, 4.5)
# </figure>

# <figure name="label_check">
wide_a <- authors_long %>% pivot_wider(names_from = sample_id, values_from = authors, values_fill = 0)
wide_o <- ours_long %>% pivot_wider(names_from = sample_id, values_from = ours, values_fill = 0)
shared <- intersect(wide_a$genus, wide_o$genus)
A <- as.matrix(wide_a[match(shared, wide_a$genus), -1])
O <- as.matrix(wide_o[match(shared, wide_o$genus), -1])
cross <- cor(O, A, method = "spearman")
ord_o <- c(intersect(sample_order, rownames(cross)), setdiff(rownames(cross), sample_order))
ord_a <- c(intersect(sample_order, colnames(cross)), setdiff(colnames(cross), sample_order))
p_label <- as.data.frame(cross) %>% rownames_to_column("ours") %>%
  pivot_longer(-ours, names_to = "authors", values_to = "rho") %>%
  mutate(ours = factor(ours, levels = rev(ord_o)), authors = factor(authors, levels = ord_a)) %>%
  ggplot(aes(authors, ours, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#f4f3ef", high = "#104281", name = "Spearman rho") +
  labs(x = "Authors' table", y = "This analysis (labels from SDRF / file names)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = if (n_samples > 40) 6 else NULL),
        axis.text.y = element_text(size = if (n_samples > 40) 6 else NULL),
        panel.grid = element_blank(), legend.position = "right")
save_plot(p_label, "label_check", max(10, 3 + 0.12 * n_samples), max(6.5, 1.5 + 0.11 * n_samples))
# </figure>
best_match <- tibble(sample = rownames(cross),
                     best_authors_sample = colnames(cross)[apply(cross, 1, which.max)],
                     rho_best = round(apply(cross, 1, max), 2),
                     rho_same_label = round(vapply(rownames(cross), function(s)
                       if (s %in% colnames(cross)) cross[s, s] else NA_real_, 1), 2))
save_table(best_match, "authors_label_check")

# ---------------------------------------------------------------- run summary and session
save_table(tibble(key = names(summary_rows), value = unlist(summary_rows)), "run_summary")
writeLines(capture.output(sessionInfo()), snakemake@output[["session"]])
