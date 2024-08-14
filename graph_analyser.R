library(igraph)

# Load the interaction data from a text file
data <- read.table("../edges.tsv", header = TRUE, sep = "\t")

# Simplify the graph to remove multiple edges and self-loops, and convert it to a data frame
semRedRede <- as_data_frame(simplify(graph_from_data_frame(data, directed=FALSE)))

# Create an undirected graph from the simplified data frame
rede <- graph_from_data_frame(semRedRede, directed = FALSE)

# Calculate the betweenness centrality for each vertex in the network
gargalo <- betweenness(rede, v = V(rede), directed = FALSE, weights = NULL, normalized = FALSE)

# Convert the betweenness centrality to a table and then to a data frame
gargalo <- as.table(gargalo)
gargalTable <- as.data.frame(gargalo)

# Save the betweenness centrality data to a text file
write.table(x=gargalTable, file = "betweenness.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Display a summary of the betweenness centrality data
summary(gargalTable)

# Display a summary of the degree centrality data
tmp <-summary(gargalTable)

# Dividindo a string com base nos ":"
split_result <- strsplit(tmp[11], ":")[[1]]

# Removendo espaços em branco desnecessários
split_result <- trimws(split_result)

# Exibindo o vetor
cutoff_betweenness <- as.numeric(split_result[2])

cutoff_betweenness
######################### Degree Centrality ###################

# Calculate the degree centrality for each vertex in the network
grau <- degree(rede, v = V(rede), loops = TRUE, normalized = FALSE)

# Convert the degree centrality to a table and then to a data frame
grau <- as.table(grau)
degreeTable <- as.data.frame(grau)

# Save the degree centrality data to a text file
write.table(x=degreeTable, file = "degree.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Display a summary of the degree centrality data
tmp <-summary(degreeTable)

# Dividindo a string com base nos ":"
split_result <- strsplit(tmp[11], ":")[[1]]

# Removendo espaços em branco desnecessários
split_result <- trimws(split_result)

# Exibindo o vetor
cutoff_degree <- as.numeric(split_result[2])
###############################################################

centralid <- cbind(gargalTable, degreeTable$Freq)

colnames(centralid) <- c("source", "betweenness", "degree")

head(centralid)

h <- subset(centralid, betweenness < cutoff_betweenness  & degree > cutoff_degree)

hb <- subset(centralid, betweenness > cutoff_betweenness  & degree > cutoff_degree)

c <- subset(centralid, betweenness < cutoff_betweenness  & degree < cutoff_degree)

b <- subset(centralid, betweenness > cutoff_betweenness  & degree < cutoff_degree)


#### Hub proteins
write.table(x=h, file = "highest_h.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(h$source), "H"), file = "h.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#### Bottleneck proteins
write.table(x=b, file = "highest_b.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(b$source), "B"), file = "b.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


### Hub-bottlenecks proteins
write.table(x=hb, file = "highest_hb.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(hb$source), "HB"), file = "hb.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


### Common proteins
write.table(x=c, file = "common.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cbind(as.character(c$source), "C"), file = "c.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


######################## Module Detection ###################
# Run the fast greedy algorithm to identify community structure (modules) in the network
fc = cluster_fast_greedy(rede)
# Reference for the method: A Clauset, MEJ Newman, C Moore: Finding community structure in very large networks.

# Save the module memberships for each vertex
nos = as.data.frame(vertex_attr(rede))
c = as.data.frame(matrix(0, ncol = 2, nrow = length(nos[,1])))
c[,1] = nos[,1]
mod = as.data.frame(fc$membership)
c[,2] = mod[,1]

# Save the module data to a text file
write.table(x=c, file = "clusters_select.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

