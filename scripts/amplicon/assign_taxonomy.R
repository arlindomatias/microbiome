# Naive Bayesian classifier against SILVA 138.2 (to genus), then exact-match species.
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(dada2))

set.seed(snakemake@params[["seed"]])
seqtab <- readRDS(snakemake@input[["seqtab"]])
seqs <- colnames(seqtab)

# Optional chunking (params chunk_size > 0) keeps peak memory low on large ASV sets: each chunk
# (classifier + species) runs in a forked child, so its memory is returned when the child exits,
# and is saved to params checkpoint_dir, so a killed run resumes where it stopped.
chunk_size <- snakemake@params[["chunk_size"]]
if (is.null(chunk_size) || chunk_size <= 0) chunk_size <- length(seqs)
chunks <- split(seqs, ceiling(seq_along(seqs) / chunk_size))
ckpt <- snakemake@params[["checkpoint_dir"]]
if (!is.null(ckpt)) dir.create(ckpt, recursive = TRUE, showWarnings = FALSE)

classify <- function(x) {
  t <- assignTaxonomy(x, snakemake@input[["trainset"]],
                      multithread = snakemake@threads, tryRC = TRUE)
  addSpecies(t, snakemake@input[["species"]], tryRC = TRUE)
}

taxa <- do.call(rbind, lapply(seq_along(chunks), function(i) {
  f <- if (!is.null(ckpt)) file.path(ckpt, sprintf("chunk_%03d_of_%03d.rds", i, length(chunks)))
  if (!is.null(f) && file.exists(f)) {
    t <- readRDS(f)
    if (identical(rownames(t), chunks[[i]])) {
      cat(sprintf("Chunk %d/%d: from checkpoint\n", i, length(chunks)))
      return(t)
    }
  }
  cat(sprintf("Chunk %d/%d (%d ASVs)\n", i, length(chunks), length(chunks[[i]])))
  if (length(chunks) == 1) {
    t <- classify(chunks[[i]])
  } else {
    set.seed(snakemake@params[["seed"]] + i)
    job <- parallel::mcparallel(classify(chunks[[i]]))
    t <- parallel::mccollect(job)[[1]]
    if (inherits(t, "try-error") || is.null(t)) stop("chunk ", i, " failed: ", t)
  }
  if (!is.null(f)) saveRDS(t, f)
  t
}))
stopifnot(identical(rownames(taxa), seqs))

out <- data.frame(asv = paste0("ASV", seq_along(seqs)), taxa, sequence = seqs,
                  row.names = NULL, check.names = FALSE)
cat("Assigned fraction per rank:\n")
print(colMeans(!is.na(taxa)))
write.table(out, snakemake@output[["taxa"]], sep = "\t", quote = FALSE, row.names = FALSE)
