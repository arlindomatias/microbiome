# Assemble the phyloseq object (ASV counts, taxonomy, midpoint-rooted tree, sample metadata).
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(ape)
  library(phangorn)
})

seqtab <- readRDS(snakemake@input[["seqtab"]])
taxa <- read.delim(snakemake@input[["taxa"]], check.names = FALSE)
meta <- read.delim(snakemake@input[["metadata"]], colClasses = "character")

asv_ids <- paste0("ASV", seq_len(ncol(seqtab)))
stopifnot(identical(taxa$sequence, colnames(seqtab)), identical(taxa$asv, asv_ids))
seqs <- DNAStringSet(colnames(seqtab)); names(seqs) <- asv_ids
colnames(seqtab) <- asv_ids

# Samples excluded in the config (e.g. unusable raw data) are dropped here: they stay in the read
# tracking table but not in the phyloseq object.
exclude <- unlist(snakemake@params[["exclude"]])
if (length(exclude) > 0) {
  stopifnot(all(exclude %in% rownames(seqtab)))
  cat("Excluded samples:", exclude, "\n")
  seqtab <- seqtab[!rownames(seqtab) %in% exclude, , drop = FALSE]
  meta <- meta[!meta$sample %in% exclude, ]
  # ASV names and sequences must stay in step with asvs.fasta and taxonomy.tsv.
  stopifnot(all(colSums(seqtab) > 0))
}
rownames(meta) <- meta$sample
# Group levels in the order given in the config (first = reference group).
group_levels <- names(snakemake@params[["groups"]])
stopifnot(all(meta$group %in% group_levels))
meta$group <- factor(meta$group, levels = group_levels)
if ("week" %in% names(meta)) meta$week <- suppressWarnings(as.integer(meta$week))
stopifnot(setequal(rownames(seqtab), rownames(meta)))

tax_mat <- as.matrix(taxa[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
rownames(tax_mat) <- taxa$asv

tree <- midpoint(read.tree(snakemake@input[["tree"]]))

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), tax_table(tax_mat),
               sample_data(meta[rownames(seqtab), ]), phy_tree(tree), refseq(seqs))
print(ps)
saveRDS(ps, snakemake@output[["ps"]])

counts <- data.frame(asv = asv_ids, t(seqtab), taxa[, -1], check.names = FALSE)
write.table(counts, snakemake@output[["counts"]], sep = "\t", quote = FALSE, row.names = FALSE)
