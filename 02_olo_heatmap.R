# --- Libraries ---
library(ComplexHeatmap)
library(circlize)
library(seriation)


# --- Load objects from Script 1 ---
load("clustering_core.RData")  # brings in seq_list, window_diss_dist, hc, etc.

# --- Optimal Leaf Ordering (OLO) ---
seriation_result <- seriate(window_diss_dist, method = "OLO_average")
optimal_order    <- get_order(seriation_result)

hc_optimal <- hc
hc_optimal$order <- optimal_order

# --- Build nucleotide matrix (same function, minimal) ---
create_nucleotide_matrix <- function(sequences) {
  if (length(unique(nchar(sequences))) > 1)
    stop("All sequences must have the same length")
  m <- matrix(NA_character_, nrow = length(sequences), ncol = nchar(sequences[1]))
  for (i in seq_along(sequences)) m[i, ] <- strsplit(sequences[i], "", fixed = TRUE)[[1]]
  rownames(m) <- paste0("Seq", seq_along(sequences))
  colnames(m) <- paste0("Pos", seq_len(ncol(m)))
  m
}

nucleotide_matrix <- create_nucleotide_matrix(seq_list)

# --- JASPAR-style colors and numeric mapping ---
jaspar_colors <- c(
  "A" = "#109648",
  "C" = "#255C99",
  "G" = "#FBB900",
  "T" = "#EE2722"
)
col_fun <- colorRamp2(c(1,2,3,4),
                      c(jaspar_colors["A"], jaspar_colors["C"],
                        jaspar_colors["G"], jaspar_colors["T"]))

numeric_matrix <- matrix(0, nrow = nrow(nucleotide_matrix), ncol = ncol(nucleotide_matrix))
numeric_matrix[nucleotide_matrix == "A"] <- 1
numeric_matrix[nucleotide_matrix == "C"] <- 2
numeric_matrix[nucleotide_matrix == "G"] <- 3
numeric_matrix[nucleotide_matrix == "T"] <- 4

# --- Heatmap with OLO dendrogram ---
dend_optimal <- as.dendrogram(hc_optimal)

ht <- Heatmap(
  numeric_matrix,
  name = "Nucleotide",
  col = col_fun,
  cluster_rows = dend_optimal,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  na_col = "#EEEEEE",
  heatmap_legend_param = list(
    at = 1:4,
    labels = c("A","C","G","T"),
    title = NULL,
    direction = "horizontal",
    nrow = 1,
    grid_width  = grid::unit(12, "mm"),
    grid_height = grid::unit(8,  "mm"),
    labels_gp   = grid::gpar(fontsize = 10, fontface = "bold"),
    legend_gp   = grid::gpar(col = "#333333", lwd = 0.8)
  )
)

# --- Save (and draw) ---
pdf("nucleotide_heatmap_optimal_ordering.pdf", width = 10, height = 8)
draw(ht)
dev.off()
draw(ht)

message("Saved: nucleotide_heatmap_optimal_ordering.pdf")
