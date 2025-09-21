# --- Libraries (lean) ---
library(Biostrings)
library(Rcpp)
library(RcppParallel)
library(dynamicTreeCut)
library(cluster)   # for silhouette

# --- Parameters (keep your variable names) ---
fasta_file  <- "data/CONTENT_consensus_pmut_0.00_Hmean_1.3841_dHnorm_0.0000.fasta"
k           <- 10
window_size <- 10
exclusion_half_width_bases = 0  # Keeps a window

# --- Source the C++ dissimilarity code ---
Rcpp::sourceCpp("exclusion.cpp")  # provides calculate_window_dissimilarity_revised_cpp_parallel()

# --- Read sequences ---
dna_sequences <- readDNAStringSet(fasta_file)
seq_list <- as.character(dna_sequences)
n <- length(seq_list)

# For ARI calculation of the sample data
# true_labels <- rep(1:5, each = 2000)

# --- Compute dissimilarities (C++) ---
seq_list_for_cpp <- as.list(seq_list)

result_cpp <- calculate_window_dissimilarity_revised_cpp_parallel(
  seq_list_r = seq_list_for_cpp,
  k           = k,
  window_size = window_size,
  exclusion_half_width_bases = exclusion_half_width_bases
)

# Primary dissimilarity matrix and dist object
window_diss_matrix <- result_cpp[["window_dissimilarity_matrix"]]
basic_diss_matrix  <- result_cpp[["basic_dissimilarity_matrix"]]
window_diss_dist  <- as.dist(window_diss_matrix)

# --- Hierarchical clustering ---
hc <- hclust(window_diss_dist, method = "average")

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
    if (length(keep) < 3L || k < 2L)  # not enough for silhouette
      return(c(ms = ms, k = k, mean_sil = NA_real_))

    # Each cluster needs >=2 members
    if (any(table(lab[keep]) < 2L))
      return(c(ms = ms, k = k, mean_sil = NA_real_))

    # Subset distance for assigned cells only
    Dk <- as.dist(as.matrix(window_diss_dist)[
      keep, keep, drop = FALSE
    ])
    sil <- cluster::silhouette(lab[keep], Dk)
    mean_sil <- mean(sil[, "sil_width"], na.rm = TRUE)
    c(ms = ms, k = k, mean_sil = mean_sil)
  }

  scores_mat <- t(sapply(sizes, eval_ms))
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
    # Fallback: global maximum (treat NA as -Inf)
    pick <- which.max(replace(y, is.na(y), -Inf))
  }

  list(
    best_minClusterSize = scores$minClusterSize[pick],
    trace = scores
  )
}

res <- auto_minClusterSize()

# --- Final clustering via cutreeHybrid  ---
clusters <- cutreeHybrid(
  dendro = hc,
  distM  = window_diss_matrix,
  minClusterSize = res$best_minClusterSize,
  deepSplit = 1,
  pamRespectsDendro = FALSE
)

cluster_labels <- clusters$labels
cat("Number of clusters identified:", max(cluster_labels, na.rm = TRUE), "\n")
print(table(cluster_labels))

# --- Save for (OLO + heatmap) ---
save(list = c(
  "fasta_file","k","window_size",
  "seq_list","basic_diss_matrix","window_diss_dist","hc",
  "clusters","cluster_labels","res"
  # ,"true_labels"  # keep if you use it elsewhere
), file = "clustering_core.RData")

message("Saved: clustering_core.RData")
