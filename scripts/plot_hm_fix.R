#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-g", "--gene"), type="character", default=NULL, 
                help="best trees file (only fourth column)", dest = "gene_trees"),
    make_option(c("-i", "--info"), type="character", default=NULL, 
                help="annotated info file", dest = "meta"),
    make_option(c("-s", "--sptree"), type = "character", default=NA,
                help="rooted species tree", dest = "sptree"),
    make_option(c("-a", "--absence"), type = "logical", default=FALSE, action = "store_true",
                help="if active only presence/absence plot", dest = "absence"),
    # make_option(c("-t", "--tax"), type = "logical", default=FALSE, action = "store_true",
    #             help="if active clustering will be based on taxonomy", dest = "tax"),
    # make_option(c("-p", "--pca"), type = "logical", default=FALSE, action = "store_true",
    #             help="if active also plot pca", dest = "pca"),
    make_option(c("--subsample"), type = "integer", default=0,
                help="sample x random trees", dest = "ntrees"),
    make_option(c("-m", "--max_missing"), type="integer", default=0.7,
                help="proportion of maximum missing samples", dest = "maxmiss"),
    make_option(c("-o", "--output"), type="character", default=NA,
                help="output heatmap name", metavar="character", dest = "outfile"),
    make_option(c("--occupancy"), type="character", default=NA,
                help="output occupancy name", metavar="character", dest = "occupancy"),
    make_option(c("-e", "--exclude"), type="character", default=NA,
                help="species to exclude, i.e. outgroups", metavar="character", dest = "exclude"),
    make_option(c("-n", "--numclusters"), type="integer", default=5,
                help="number of clusters to split columns", dest = "nclust")
)

library(ape)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

col_palette <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                          "#D55E00","#CC79A7","#999999","#000000")

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

names_info <- c("taxid", "mnemo", "date", "longest", "source", "sp_og",
                "k","p","c","o","f","g","s")

info <- read.delim(opt$meta, col.names = names_info, header = F) %>% 
  mutate(sp_og=gsub("\\]", "", gsub("\\[", "", sp_og)), label = sub("^([A-Za-z])[A-Za-z]+\\s+([A-Za-z]+)", "\\1. \\2", sp_og))
info$mnemo <- gsub("\\..*", "", info$mnemo)

# read gene_trees
ptrees <- read.tree(opt$gene_trees)
if (opt$ntrees != 0){
    ptrees <- sample(ptrees, opt$ntrees)
}

if (!is.na(opt$exclude)){
  remove <- unlist(strsplit(opt$exclude, " "))
  } else {
  remove = c()
}

# read species tree and transform in dendrogram
# if not sptree given then cluster by hclust
if (!is.na(opt$sptree)){
    species_tree <- ladderize(drop.tip(read.tree(opt$sptree), remove))
    species <- unique(species_tree$tip.label)
    if (is.null(species_tree$edge.length)) {
      species_tree$edge.length <- rep(1, length(species_tree$node.label) + length(species_tree$tip.label) - 1)
    }
    row_clust <- as.hclust(chronos(species_tree))
} else {
    species <- unique(unlist(sapply(ptrees, function(x) gsub(".*_" , "", x$tip.label), simplify = T)))
}

n_cat <- 9
unique_tax <- abs(as.numeric(apply(info[, 7:12], 2, function(x) length(unique(x)))) - n_cat)
names(unique_tax) <- colnames(info)[7:12]
grp_col <- names(which(unique_tax==min(unique_tax)))[1]

add_colors = grDevices::colors()[grep('gr(a|e)y|black', grDevices::colors(), invert = T)]
miss_colors <- length(unique(info[, grp_col])) - length(col_palette)

if (miss_colors>0) {
  col_palette <- c(col_palette, add_colors[1:miss_colors])
}

names(col_palette) <- unique(info[, grp_col])
col_palette <- col_palette[which(!is.na(names(col_palette)))]
# print(col_palette)
species <- setdiff(species, remove)

min_sp <- floor(opt$maxmiss*length(species))

matrix_presence <- do.call("rbind", (lapply(ptrees, function(x) sapply(species, function(sp) sum(grepl(sp, x$tip.label))))))

# matrix_presence[which(matrix_presence==0)] <- NA

max_abs = quantile(matrix_presence, .95)

max_col = ifelse(max(matrix_presence, na.rm = T)>max_abs, max_abs, max(matrix_presence, na.rm = T))
steps = max_col/4


if (opt$absence) {
    col_fun = circlize::colorRamp2(c(0, 1), c("white", "black"))
    matrix_presence <- matrix_presence[apply(matrix_presence, 1, function(x) max(x)<2), ]
    title = "Present"
    # matrix_presence[which(matrix_presence>0)] <- 1
} else {
    breaks <- c(0, 1, seq(2, max_col, steps))
    cols <- c("white", "black",  c(viridis::plasma(n = 4)))[1:length(breaks)]
    col_fun = circlize::colorRamp2(breaks, cols)#as.character(wesanderson::wes_palette("Zissou1"))))
    title = "Size"
}

if (is.na(opt$sptree)) {
    row_clust <- dendsort::dendsort(fastcluster::hclust(dist(t(matrix_presence), method = "man")))
}

grp <- info[match(colnames(matrix_presence), info$mnemo), grp_col]

width = log10(dim(matrix_presence)[1])*3.5
height = dim(matrix_presence)[2]/3.8

# genes_cluster = dendextend::color_branches(genes_cluster, k = 5)
fh = function(x) fastcluster::hclust(dist(x, method = "man"))

# outmat <- ifelse(opt$absence, paste0(opt$outfile, "_presence.png"), paste0(opt$outfile, "_spcounts.png"))
outmat <- opt$outfile
pass <- rowSums(matrix_presence == 1) > min_sp & rowSums(matrix_presence > 1) == 0

pdf(outmat, width = 25, height = 8)

ha = HeatmapAnnotation(
  pass = pass, col = list(pass=c("TRUE"="black", "FALSE"="white")),
  annotation_name_side = "left", 
  annotation_legend_param = list(pass = list(direction="horizontal")))

hr = rowAnnotation(clade=grp, col = list(clade=col_palette), annotation_name_side = "top",
                   annotation_legend_param = list(clade = list(nrow=3, direction="horizontal")))

new_labels <- info[match(colnames(matrix_presence), info$mnemo), "s"]


hm = Heatmap(t(matrix_presence),
            title,
            cluster_columns = fh,
            column_km = opt$nclust,
            use_raster = TRUE,
            row_names_max_width = unit(150, "cm"),
            raster_quality = 5,
            column_title = NULL,
            column_dend_gp = gpar(lwd = .2),
            rect_gp = gpar(col = NA, lwd = 0),
            border = TRUE,
            cluster_rows = row_clust,
            col = col_fun, 
            top_annotation = ha,
            right_annotation = hr,
            heatmap_legend_param = list(direction = "horizontal"))

ht_list = hm + 
  rowAnnotation(labels = anno_text(new_labels, which = "row"))


draw(ht_list, merge_legend = TRUE,  heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


# if (opt$pca){
#   print("doing pca")
#     pca <- prcomp(t(matrix_presence), center = TRUE)
#     out_pca <- paste0(opt$outfile,"_pca.png")
#     plot_pca <- factoextra::fviz_pca_ind(pca, )
#    ggplot2::ggsave(out_pca, plot_pca, device = "png", width = 7, height = 7, dpi=300)
# } else{
#   print("done")
# }


if (!is.na(opt$occupancy)) {
  n_genes <- nrow(matrix_presence)
  n_sp <- ncol(matrix_presence)
  
  # theme_set(theme_bw(base_family = "Helvetica"))
  theme_set(theme_bw())
  
  df_sp <- tibble(mnemo=colnames(matrix_presence),
                  occupancy=colSums(matrix_presence!=0)/n_genes,
                  expanded=colSums(matrix_presence>1)/n_genes)
  
  occupancy <- left_join(df_sp, info) %>% 
    ggplot(aes(occupancy, expanded)) + 
    geom_point(aes(fill=.data[[grp_col]]), pch=21, alpha=0.9, size=2) + 
    scale_fill_manual(values = col_palette) + 
    ggrepel::geom_text_repel(aes(label=mnemo), size=2) +
    labs(x="occupancy", y="expanded", fill="clade")
  
  
  ggsave(opt$occupancy, occupancy, width = 8, height = 7)
}
