suppressPackageStartupMessages(library(ComplexHeatmap))
suppressMessages(library(optparse))
suppressMessages(library(igraph))

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="diagonal matrix", dest = "input"),
  make_option(c("-m","--minsize"), type = "integer", default=4,
              help="minimum cluster size. Default 4", dest = "minclust"),
  make_option(c("-p","--perc"), type = "double", default=.25,
              help="minimum proportion of shared homologs. Default .25", dest = "perc"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output plot", dest = "output")
  # make_option(c("-a","--approximate"), type = "logical", default=FALSE, action = "store_true",
  #             help="if active, reduce matrix by approximate method", dest = "approx")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(is.null(opt$input) | is.null(opt$output)) { 
  cat("Both -i and -o are required!\n")
  print_help(opt_parser)
  quit(status=1)
}

# opt <- NULL
# opt$input <- "a_distance.tsv"
# opt$output <- "plots/a.pdf"
# opt$minclust <- 50
# opt$perc <- .8
# pass <- readLines("outputs/ids.txt")


df_mat <- data.table::fread(opt$input, data.table = FALSE)
rownames(df_mat) <- df_mat[,1]
df_mat <- as.matrix(df_mat[, -1])
n_homo <- diag(df_mat)

# df_mat <- df_mat[pass, pass]

# remove those with no overlapping sets
max_0 <- nrow(df_mat)-1
idxs <- which(rowSums(df_mat==0)<max_0)

idxs_uniq <- which(rowSums(df_mat==0)==max_0)

print(summary(rowSums(df_mat!=0)))
print(paste(length(idxs_uniq), "out of", nrow(df_mat), "are unique"))
# print(length(which(rowSums(df_mat!=0)>opt$minclust)))

df_mat <- df_mat[names(idxs), names(idxs)]

print("scaling matrix...")
make_symmetric <- function(m, upper=TRUE) {
  out_m <- m
  diag_m <- out_m*diag(nrow(out_m))
  if (upper) {
    out_m[upper.tri(out_m, diag = TRUE)] <- 0
  } else {
    out_m[lower.tri(out_m, diag = TRUE)] <- 0
  }
  out_m <- out_m + t(out_m) + diag_m
  # out_m + t(out_m)
  # out_m[upper.tri(out_m)] <- t(out_m)[upper.tri(out_m)]
  return(out_m)
}

df_mat <- t(t(df_mat) / diag(df_mat))
df_mat_lo <- make_symmetric(df_mat, upper = FALSE)
df_mat_up <- make_symmetric(df_mat)
df_mat <- purrr::reduce(list(df_mat_lo,df_mat_up),`+`) / 2
# df_mat <- apply(df_mat, 1, scales::rescale)
df_mat <- df_mat/diag(df_mat)
print("scaling done!")

# easiest way I found to compute disjoint sets
print("filtering sets...")

df_mat_adj <- df_mat
# if you remove low percentage identities from dataset
# subnetworks are not easily influenced by a single hit!

df_mat_adj[df_mat_adj<opt$perc] <- 0
df_mat_adj[t(df_mat_adj<opt$perc)] <- 0

adj_network <- graph_from_adjacency_matrix(df_mat_adj, weighted = T)
split_graph <- decompose.graph(adj_network, mode = "strong", min.vertices = opt$minclust)
n_subgraphs <- length(split_graph)
print(n_subgraphs)
lengths_graphs <- sapply(split_graph, length)
summary(lengths_graphs)
idxs_good <- V(disjoint_union(split_graph))$name
print(paste("filtering done:", length(idxs_good), "will be analyzed"))


df_mat <- df_mat[idxs_good, idxs_good]

# if (opt$approx){
#   print("approximate method: some seeds may lose corresponding hits!")
#   idxs_good <- which(rowSums(df_mat!=0)>opt$minclust)
# 
#   } else {
#   print("first clustering...")
#   clust_mat <- fastcluster::hclust(dist(df_mat, "binary"), "ave")
#   print("first clustering done!")
#   
#   cut_clust <- cutree(clust_mat, h=.99999) 
#   idxs_good <- names(cut_clust[cut_clust %in% which(table(cut_clust) > opt$minclust)])
#   }

print(paste("Computing matrix for:", nrow(df_mat)))

print("final clustering...")
clust_mat <- fastcluster::hclust(dist(df_mat), "complete")
print("final clustering done!")

# ifelse(opt$perc>.51, .51, opt$perc)

if (opt$perc<=1){
  colors <- c("black", "#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00", "white")
  breaks <- c(0, seq(0.1, 0.99, length.out=5), 1)
} else {
  colors <- c("black", "grey", "#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00", "white")
  breaks <- c(0, opt$perc, seq(opt$perc+1e-6, .99, length.out=5), 1)
}

col_fun = circlize::colorRamp2(breaks, colors)

n_homo <- n_homo[idxs_good]
column_ha = HeatmapAnnotation(n_homo = anno_barplot(n_homo))

ht = Heatmap(df_mat, name = "shared homologs",
             col = col_fun, 
             cluster_columns = clust_mat,
             cluster_rows = clust_mat,
             row_split = ifelse(n_subgraphs==1, 2, n_subgraphs), 
             row_title = NULL,
             column_split = ifelse(n_subgraphs==1, 2, n_subgraphs),
             column_title = NULL,
             row_gap = unit(c(0), "mm"),
             column_gap = unit(c(0), "mm"),
             # row_gap = unit(c(2, 4), "mm"),
             show_column_dend = F, show_row_dend = F,
             top_annotation = column_ha,
             rect_gp = gpar(col = NA, lwd = 0),
             # show_heatmap_legend = FALSE,
             row_labels = rep("", nrow(df_mat)),
             column_labels = rep("", nrow(df_mat)),
             # use_raster = FALSE,
             use_raster = TRUE, raster_quality = 5,
             heatmap_legend_param = list(
             at = c(0, 0.5, 0.9, 1),
             legend_direction = "horizontal", 
             legend_width = unit(8, "cm"), 
             legend_heigth = unit(.5, "cm")))
# draw(ht)
# w = ComplexHeatmap:::width(ht)
# w = convertX(w, "inch", valueOnly = TRUE)
# h = ComplexHeatmap:::height(ht)
# h = convertY(h, "inch", valueOnly = TRUE)
dims = nrow(df_mat)*1.5/100
dims = ifelse(dims < 10, 10, dims)
dims = ifelse(dims > 20, 20, dims)

pdf(file=opt$output, width = dims, height = dims)
draw(ht, heatmap_legend_side="bottom")
dev.off()

# longest_graph <- which(lengths_graphs==max(lengths_graphs))[1]
# writeLines(V(split_graph[[longest_graph]])$name, "outputs/ids.txt")
