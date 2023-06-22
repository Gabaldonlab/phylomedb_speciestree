library('TreeDist')

opt <- NULL
opt$gene_trees <- "output/data/all_trees/best_trees_0003_labs.nwk"
opt$ntrees <- 100
opt$maxmiss <- .7
opt$meta <- "output/data/info/info_0003.tsv"
opt$outfile <- "test.png"


names_info <- c("taxid", "mnemo", "date", "longest", "source", "sp_og",
                "k","p","c","o","f","g","s")

info <- read.delim(opt$meta, col.names = names_info, header = F)
info$mnemo <- gsub("\\..*", "", info$mnemo)

# read gene_trees
ptrees <- read.tree(opt$gene_trees)

matrix_presence <- do.call("rbind", (lapply(ptrees, function(x) sapply(species, function(sp) sum(grepl(sp, x$tip.label))))))

min_sp <- floor(opt$maxmiss*length(species))

pass <- rowSums(matrix_presence == 1) > min_sp & rowSums(matrix_presence > 1) == 0

ptrees <- ptrees[pass]

distances <- ClusteringInfoDistance(ptrees)
sammon <- MASS::sammon(distances, k = 12)
mapping <- sammon$points
plot(mapping,
     asp = 1, # Preserve aspect ratio - do not distort distances
     ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
)
