#!/usr/bin/env Rscript

# --- Libraries (silent) ---
suppressPackageStartupMessages({
  library(Biostrings)
  library(Rcpp)
  library(RcppParallel)
  library(dynamicTreeCut)
  library(cluster)   
})

# --- CLI parameters (mandatory fasta, k, window; optional excl=0) ---
args <- commandArgs(trailingOnly = TRUE)

# tiny --key=val parser
kv <- if (length(args)) {
  parts <- strsplit(args, "=", fixed = TRUE)
  setNames(vapply(parts, `[`, "", 2L), vapply(parts, `[`, "", 1L))
} else character()

# usage helper
usage <- function() {
  cat("Usage:\n",
      "  Rscript clustering.R --fasta=FILE --k=INT --window=INT [--excl=INT]\n",
      "  Required: --fasta, --k, --window\n",
      "  Optional: --excl (default 0)\n", sep = "")
}

# helper
if (!is.na(kv["--help"]) || any(args %in% c("-h", "--help"))) {
  usage(); quit(status = 0)
}

# validate presence
if (is.na(kv["--fasta"]))  stop("Missing required argument: --fasta=<file>\n", call. = FALSE)
if (is.na(kv["--k"]))      stop("Missing required argument: --k=<integer>\n", call. = FALSE)
if (is.na(kv["--window"])) stop("Missing required argument: --window=<integer>\n", call. = FALSE)

# assign + validate type
fasta_file <- kv["--fasta"]

k <- as.integer(kv["--k"])
if (is.na(k) || k < 1) stop("`--k` must be a positive integer\n", call. = FALSE)

window_size <- as.integer(kv["--window"])
if (is.na(window_size) || window_size < 1) stop("`--window` must be a positive integer\n", call. = FALSE)

exclusion_half_width_bases <- if (!is.na(kv["--excl"])) as.integer(kv["--excl"]) else 0
if (is.na(exclusion_half_width_bases) || exclusion_half_width_bases < 0)
  stop("`--excl` must be a non-negative integer (default 0)\n", call. = FALSE)

# --- Mini logger (timestamps) ---
stamp <- function(fmt, ...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      sprintf(fmt, ...), "\n", sep = "")
}

stamp("Params => fasta=%s | k=%d | window=%d | excl=%d",
      fasta_file, k, window_size, exclusion_half_width_bases)

# --- Sourcing the C++ code for computing dissimilarity matrix  ---
suppressMessages(Rcpp::sourceCpp("exclusion.cpp"))  # provides calculate_window_dissimilarity_revised_cpp_parallel()

# --- Read sequences ---
stamp("Reading FASTA: %s", fasta_file)
dna_sequences <- readDNAStringSet(fasta_file)
seq_list <- as.character(dna_sequences)
n <- length(seq_list)
stamp("Loaded %d sequences", n)

# --- Compute dissimilarities (C++) ---
stamp("Computing dissimilarities in C++ (k=%d, window=%d, excl=%d) ...",
      k, window_size, exclusion_half_width_bases)
t0 <- proc.time()
seq_list_for_cpp <- as.list(seq_list)

result_cpp <- calculate_window_dissimilarity_revised_cpp_parallel(
  seq_list_r = seq_list_for_cpp,
  k           = k,
  window_size = window_size,
  exclusion_half_width_bases = exclusion_half_width_bases
)
stamp("Dissimilarities done in %.1f sec", (proc.time() - t0)[["elapsed"]])

# Primary dissimilarity matrix and dist object
window_diss_matrix <- result_cpp[["window_dissimilarity_matrix"]]
basic_diss_matrix  <- result_cpp[["basic_dissimilarity_matrix"]]
stamp("Building dist object ...")
window_diss_dist <- as.dist(window_diss_matrix)

# --- Hierarchical clustering ---
stamp("Running hclust (average) ...")
t_hc <- proc.time()
hc <- hclust(window_diss_dist, method = "average")
stamp("hclust done in %.1f sec", (proc.time() - t_hc)[["elapsed"]])

# --- Minimal auto_minClusterSize: mean silhouette plateau; else global max ---
auto_minClusterSize <- function(step = 0.025, max_frac = 0.5,
                                deepSplit = 1, w = 3, eps = 0.003,
                                pamStage = FALSE) {
  stopifnot(inherits(window_diss_dist, "dist"))
  N <- attr(window_diss_dist, "Size")

  # Candidate grid (at least 2)
  sizes <- sort(unique(pmax(2L, round(seq(step, max_frac, by = step) * N))))
  if (!length(sizes)) sizes <- 2L

  eval_ms <- function(ms) {
    lab <- try(
      dynamicTreeCut::cutreeHybrid(
        dendro = hc, distM = window_diss_matrix,
        minClusterSize = ms, deepSplit = deepSplit,
        pamStage = pamStage, pamRespectsDendro = FALSE
      )$labels,
      silent = TRUE
    )
    if (inherits(lab, "try-error") || is.null(lab))
      return(c(ms = ms, k = NA_real_, mean_sil = NA_real_))

    keep <- which(lab > 0L)
    k <- if (length(keep)) length(unique(lab[keep])) else 0L
    if (length(keep) < 3L || k < 2L)
      return(c(ms = ms, k = k, mean_sil = NA_real_))

    if (any(table(lab[keep]) < 2L))
      return(c(ms = ms, k = k, mean_sil = NA_real_))

    Dk <- as.dist(as.matrix(window_diss_dist)[keep, keep, drop = FALSE])
    sil <- cluster::silhouette(lab[keep], Dk)
    mean_sil <- mean(sil[, "sil_width"], na.rm = TRUE)
    c(ms = ms, k = k, mean_sil = mean_sil)
  }

  stamp("Scanning minClusterSize grid (%d candidates) ...", length(sizes))
  pb <- txtProgressBar(min = 1, max = length(sizes), style = 3)
  scores_mat <- matrix(NA_real_, nrow = length(sizes), ncol = 3)
  for (i in seq_along(sizes)) {
    scores_mat[i, ] <- eval_ms(sizes[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)

  colnames(scores_mat) <- c("ms", "k", "mean_sil")
  scores <- data.frame(
    minClusterSize   = as.integer(scores_mat[, "ms"]),
    k                = as.integer(scores_mat[, "k"]),
    mean_silhouette  = as.numeric(scores_mat[, "mean_sil"]),
    row.names = NULL
  )

  # Plateau detection on mean silhouette
  y  <- scores$mean_silhouette
  dy <- c(NA_real_, diff(y))
  stable <- vapply(seq_along(y), function(i) {
    if (i <= w) return(FALSE)
    wnd <- dy[(i - w + 1):i]
    if (any(is.na(wnd))) return(FALSE)
    max(abs(wnd)) < eps
  }, logical(1))

  pick <- which(stable & scores$k >= 2L)[1]
  if (is.na(pick)) {
    pick <- which.max(replace(y, is.na(y), -Inf))
  }

  stamp("Chosen minClusterSize = %d", scores$minClusterSize[pick])
  list(
    best_minClusterSize = scores$minClusterSize[pick],
    trace = scores
  )
}

# --- Setting the minClusterSize to auto ---
res <- auto_minClusterSize()

# --- Final clustering via cutreeHybrid ---
stamp("Cutting tree with cutreeHybrid ...")
clusters <- cutreeHybrid(
  dendro = hc,
  distM  = window_diss_matrix,
  minClusterSize = res$best_minClusterSize,
  deepSplit = 1,
  pamRespectsDendro = FALSE
)
cluster_labels <- clusters$labels
stamp("Identified %d clusters", max(cluster_labels, na.rm = TRUE))

# Optional: quick cluster size table (concise)
tbl <- table(cluster_labels)
print(tbl)

# --- Save for (OLO + heatmap) ---
stamp("Saving clustering_core.RData ...")
save(list = c(
  "fasta_file","k","window_size",
  "seq_list","basic_diss_matrix","window_diss_dist","hc",
  "clusters","cluster_labels","res"
), file = "clustering_core.RData")
stamp("Saved: clustering_core.RData")
