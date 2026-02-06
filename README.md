# Network Analysis: Hubs, Bottlenecks, and Community Detection

This repository contains an R-based pipeline for Protein-Protein Interaction (PPI) network analysis. It calculates topological metrics to identify key proteins (Hubs and Bottlenecks) and detects functional modules within the network.

## Description

The script processes interaction data to:
1.  **Calculate Centrality:** Computes Betweenness and Degree for each node.
2.  **Classify Proteins:** Categorizes nodes based on automated thresholds (usually the Mean or Median) into:
    * **Hubs (H):** High degree, low betweenness (highly connected nodes).
    * **Bottlenecks (B):** Low degree, high betweenness (nodes that bridge different clusters).
    * **Hub-Bottlenecks (HB):** High degree and high betweenness (key regulatory nodes).
    * **Common (C):** Low degree and low betweenness.
3.  **Community Detection:** Identifies network modules using the Fast-Greedy algorithm.
4.  **Data Organization:** Automatically exports all results into a structured folder.



## Dependencies

### Software Requirements
* **R (>= 4.0.0)**

### R Packages
The core analysis is powered by the `igraph` package. You can install it via R console:

```r
install.packages("igraph")
```

### Project Structure

`analysis.R`: The main script containing the pipeline.

`exemple/`: Directory for input files (e.g., STRING database exports).

`results/`: Directory automatically created by the script to store outputs.

### Usage

Prepare Data: Place your interaction file in the exemple/ folder. The default script expects a space-separated file, but you can adjust the sep parameter in read.table.

```r
data <- read.table("exemple/your_data.txt", header = TRUE, sep = " ")
```

### Output

All results are saved in the `results/` folder:

`betweenness.txt`: Calculated betweenness centrality for all nodes.

`degree.txt`: Calculated degree centrality for all nodes.

`highest_h.txt`: List of identified Hub proteins.

`highest_b.txt`: List of identified Bottleneck proteins.

`highest_hb.txt`: List of identified Hub-Bottleneck proteins.

`clusters_select.txt`: Mapping of nodes to their respective modules/clusters.

`cluster_N.txt`: Individual list of proteins for each detected module.
