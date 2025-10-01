# InstaClust

Code and Data used for the InstaClust framework for R

## What it does

InstaClust is a fast and deterministic clustering algorithm designed to identify regulatory DNA motifs by grouping related sequences using a simple Hamming distance and hierarchical clustering. 

1.  **Compute Clusters**: The first script, `01_compute_clusters.R`, takes a FASTA file containing DNA sequences and performs hierarchical clustering on them. It calculates a dissimilarity matrix based on k-mer distances and then uses the `cutreeHybrid` function to identify clusters. The results, including the clustering information and dissimilarity matrix, are saved to a file named `clustering_core.RData`.

2.  **Generate Heatmap**: The second script, `02_olo_heatmap.R`, loads the `clustering_core.RData` file. It then uses Optimal Leaf Ordering (OLO) to reorder the sequences within the clustering for visualization. Finally, it generates a heatmap of the nucleotide sequences.
## Required Packages

To run InstaClust, you will need the following R packages installed:

* `Biostrings`
* `Rcpp`
* `RcppParallel`
* `dynamicTreeCut`
* `cluster`
* `ComplexHeatmap`
* `circlize`
* `seriation`

## How to Run
InstaClust is run from the command line in two steps:

### Step 1: Compute Clusters

Run the `01_compute_clusters.R` script with your FASTA file and parameters:

```bash
Rscript 01_compute_clusters.R --fasta=data/CONTENT_consensus_pmut_0.00_Hmean_1.3841_dHnorm_0.0000.fasta --k=10 --window=12
```

Arguments 

--fasta: Path to your input FASTA file.

--k: The k-mer size to use for distance calculation.

--window: The window size for the dissimilarity calculation.

--excl (optional): The half-width of the central exclusion region in bases. Defaults to 0.

### Step 2: Generate Heatmap

After the first step is complete and clustering_core.RData has been created, run the 02_olo_heatmap.R script.

```bash
Rscript 02_olo_heatmap.R
```
This script will load the clustering_core.RData file and generate the heatmap file in the same directory:
