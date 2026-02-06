library(igraph)

# Create the results directory if it doesn't already exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Load the interaction data from a text file
# Note: Ensure the path to the input file is correct for your environment
data <- read.table("exemple/STRG0A37MEC.protein.links.v12.0_400.txt", header = TRUE, sep = " ")

# Simplify the graph to remove multiple edges and self-loops, and convert it to a data frame
semRedRede <- as_data_frame(simplify(graph_from_data_frame(data, directed=FALSE)))

# Create an undirected graph from the simplified data frame
rede <- graph_from_data_frame(semRedRede, directed = FALSE)

# --- Betweenness Centrality ---

# Calculate the betweenness centrality for each vertex
gargalo <- betweenness(rede, v = V(rede), directed = FALSE, weights = NULL, normalized = FALSE)
gargalo <- as.table(gargalo)
gargalTable <- as.data.frame(gargalo)

# Save the betweenness centrality data to the results folder
write.table(x=gargalTable, file = "results/betweenness.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Extract the cutoff value (Mean/Median) from the summary string logic
tmp <- summary(gargalTable)
split_result <- trimws(strsplit(tmp[11], ":")[[1]])
cutoff_betweenness <- as.numeric(split_result[2])

# --- Degree Centrality ---

# Calculate the degree centrality for each vertex
grau <- degree(rede, v = V(rede), loops = TRUE, normalized = FALSE)
grau <- as.table(grau)
degreeTable <- as.data.frame(grau)

# Save the degree centrality data to the results folder
write.table(x=degreeTable, file = "results/degree.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Extract the cutoff value for degree
tmp2 <- summary(degreeTable)
split_result_deg <- trimws(strsplit(tmp2[11], ":")[[1]])
cutoff_degree <- as.numeric(split_result_deg[2])

# --- Protein Classification ---

centralid <- cbind(gargalTable, degreeTable$Freq)
colnames(centralid) <- c("source", "betweenness", "degree")

# Classify proteins based on thresholds
h  <- subset(centralid, betweenness < cutoff_betweenness & degree > cutoff_degree)
hb <- subset(centralid, betweenness > cutoff_betweenness & degree > cutoff_degree)
c_comm  <- subset(centralid, betweenness < cutoff_betweenness & degree < cutoff_degree)
b  <- subset(centralid, betweenness > cutoff_betweenness & degree < cutoff_degree)

# Save classified proteins to the results folder
# Hubs
write.table(x=h, file = "results/highest_h.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(h$source), "H"), file = "results/h.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Bottlenecks
write.table(x=b, file = "results/highest_b.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(b$source), "B"), file = "results/b.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Hub-Bottlenecks
write.table(x=hb, file = "results/highest_hb.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(hb$source), "HB"), file = "results/hb.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Common
write.table(x=c_comm, file = "results/common.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(c_comm$source), "C"), file = "results/c.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# --- Module Detection ---

# Run the fast greedy algorithm to identify community structure
fc = cluster_fast_greedy(rede)
nos = as.data.frame(vertex_attr(rede))
cluster_map = as.data.frame(matrix(0, ncol = 2, nrow = length(nos[,1])))
cluster_map[,1] = nos[,1]
cluster_map[,2] = as.data.frame(fc$membership)[,1]

# Save the module data to the results folder
write.table(x=cluster_map, file = "results/clusters_select.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# --- Individual Cluster Export ---

unique_clusters <- unique(cluster_map[,2])

for (cluster in unique_clusters) {
  # Filter nodes belonging to the current cluster
  subset_cluster <- cluster_map[cluster_map[,2] == cluster, 1, drop = FALSE]
  
  # Define file path inside the results folder
  file_name <- paste0("results/cluster_", cluster, ".txt")
  
  # Save each cluster to its own file
  write.table(subset_cluster, file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}