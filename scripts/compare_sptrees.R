
suppressMessages(library(optparse))

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="sptrees separated by comma", dest = "input"),
    make_option(c("-c", "--code"), type="character", default=NULL,
                help="code phylomedb", dest = "code"),
    make_option(c("-t", "--tax"), type="character", default=NULL,
                help="pdb info file", dest = "tax"),
    make_option(c('-s', '--sptree'), type = "character", default = NULL, 
                help = "sptree to root", dest = "sptree"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="output species tree", dest = "output"),
    make_option(c("-n", "--ncat"), type="numeric", default=7, 
                help="number of taxonomic categories to divide colors", dest = "ncat")
    # make_option(c("-a","--approximate"), type = "logical", default=FALSE, action = "store_true",
    #             help="if active, reduce matrix by approximate method", dest = "approx")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressMessages(library(tidyverse))
suppressMessages(library(TreeTools))
suppressMessages(library(ggtree))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ape))
suppressMessages(library(phangorn))
# suppressMessages(library(ggtree))
# suppressMessages(library(TreeDist))

# create Results directory
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = T)

prefix <- gsub("_rooted_sptree.nwk", "", opt$output)

# Read phylome metadata to annotate sptrees
names_info <- c("taxid", "mnemo", "date", "longest", "source", "sp_og",
                "k","p","c","o","f","g","s")
info <- read.delim(opt$tax, col.names = names_info, header = F) %>% 
    mutate(sp_og=gsub("\\]", "", gsub("\\[", "", sp_og)), 
           short_s = sub("^([A-Za-z])[A-Za-z]+\\s+([A-Za-z]+)", "\\1. \\2", sp_og))
info$mnemo <- gsub("\\..*", "", info$mnemo)

col_palette <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                          "#D55E00","#CC79A7","#999999","#000000")

unique_tax <- abs(as.numeric(apply(info[, 7:12], 2, function(x) length(unique(x)))) - opt$ncat)
names(unique_tax) <- colnames(info)[7:12]
grp_col <- names(which(unique_tax==min(unique_tax)))[1]

add_colors = grDevices::colors()[grep('gr(a|e)y|black', grDevices::colors(), invert = T)]
miss_colors <- length(unique(info[, grp_col])) - length(col_palette)

if (miss_colors>0) {
  col_palette <- c(col_palette, add_colors[1:miss_colors])
}

names(col_palette) <- unique(info[, grp_col])
col_palette <- col_palette[which(!is.na(names(col_palette)))]

# Species tree files
tree_files <- unlist(strsplit(opt$input, ","))

trees <- list()
raw_trees <- list()

# Initialize root txt file for duptree and SpeciesRax
root_file <- paste0(prefix, "_roots.txt")
cat("method\toutgroup", file=root_file, append=FALSE, sep = "\n")

methods <- c()
for (file in tree_files){
    method <- gsub(".nwk", "", gsub(".*_", "", file))
    tree <- ape::read.tree(file)

    if (method %in% c("dt", "sprax")) {
        # Get root from these methods
        first_clades <- Descendants(tree, tree$edge[1], type = "children")[1:2]
        og_node <- first_clades[which.min(sapply(first_clades, function(x) length(unlist(Descendants(tree, x, type="tips")))))]
        ogs <- tree$tip.label[phangorn::Descendants(tree, og_node, type = "tips")[[1]]]

        cat(paste(trimws(method), "\t", ogs), file=root_file, append=TRUE, sep="\n")

        if (method=="dt") {
            final_og <- ogs[1]
        }
    }
    # This is to visualize all the raw species trees as outputted by the software
    if (!"edge.length" %in% names(tree)) {
        tree$edge.length <- NA
    } else {
        tree$edge.length[is.na(tree$edge.length)] <- 0
    }
    raw_trees[[method]] <- tree

    # Remove root and branch length info from trees in order to do the densitree
    tree <- unroot(tree)
    tree$edge.length <- NULL
    trees[[method]] <- tree
    methods <- c(methods, method)
    # reorder(TreeTools::RootTree(tree, og))
}

rename_tips <- function(tree, meta) {
    tree$tip.label <- meta[match(tree$tip.label, meta$mnemo), ]$short_s
    return(tree)
}

class(trees) <- 'multiPhylo'
trees <- reorder(RootTree(trees, final_og))

# RF Matrix heatmap

rfs <- RF.dist(trees, normalize = T)
rf_matrix <- as.matrix(1 - Matrix::forceSymmetric(as.matrix(rfs), uplo="L"))

rf_df <- tibble::as_tibble(rf_matrix)
rf_df$method <- rownames(rf_matrix)
rf_df <- rf_df[,c(ncol(rf_df),1:ncol(rf_df)-1)]
readr::write_delim(rf_df, paste0(prefix, "_rf.txt"), delim="\t")

ht <- Heatmap(rf_matrix, 
        "RF_methods",
        col = circlize::colorRamp2(seq(0,1,length.out=9), RColorBrewer::brewer.pal(9, "YlGnBu")),
        width = ncol(rf_matrix)*unit(.4, "in"), 
        height = nrow(rf_matrix)*unit(.4, "in"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(rf_matrix[i, j])){
                grid.text(sprintf("%1.2f", rf_matrix[i, j]), x, y, 
                          gp = gpar(fontsize = 8, col=ifelse(rf_matrix[i, j]<.75, "black", "white")))
            } else {
                NULL
            }
        })

pdf(file = paste0(prefix, "_rf.pdf"), width=5.35, height=5.35)
draw(ht)
dev.off()

# Get consensus, rooted with Duptree outgroup
consensus <- reorder(root(read.tree(opt$sptree), final_og, resolve.root = T))

if (!is.binary(consensus)) {
    paste("Consensus species tree is not binary. Multifurcations will be solved randomly")
    consensus <- multi2di(consensus)
}

# Consensus with informative species names
consensus_nms <- rename_tips(consensus, info)

# Max root2tip distance
max_x <- max(vcv(consensus))
width_plot <- ifelse(max_x*4<=8, 8, max_x*4)

class(raw_trees) <- 'multiPhylo'
names(raw_trees) <- methods

multree <- ggtree(raw_trees) + 
    facet_wrap(~.id, scales = "free_x", ncol = 4) + 
    geom_tiplab(size=2) +
    geom_nodelab(aes(x=branch, label=round(as.numeric(label), 2)), 
                 nudge_y = .3, color="red", size=1) +
    geom_treescale() +
    scale_x_continuous(expand = c(.2,.2))

ggsave(paste0(prefix, "_multree.pdf"), multree, units="in", width=12, height=8)


trees <- sapply(trees, function(x) rename_tips(x, info))
class(trees) <- 'multiPhylo'
names(trees) <- methods

n_seqs_coeff <- length(consensus$tip.label)/4
height_plot <- ifelse(n_seqs_coeff<=5, 5, n_seqs_coeff)
# offset_text <- length(consensus$tip.label)/80

df_cons <- tibble::as_tibble(consensus) %>% 
    left_join(info, by=c("label"="mnemo")) %>% 
    separate(label, c("supname", "gCF", "sCF"), sep = "/") %>% 
    rowwise() %>% 
    mutate(support=ifelse(as.numeric(supname)<=100, as.numeric(supname), NA),
           gCF=as.numeric(gCF),
           sCF=as.numeric(sCF),
           mean_sup = mean(c(support, gCF, sCF), na.rm=T),
           tiplab=ifelse(is.na(support), short_s, ""))

right_x <- max_x+max(nchar(consensus_nms$tip.label))*(0.02*max_x)

cons_plot <- ggtree(consensus) %<+% df_cons +
    geom_tiplab(aes(label=paste0(" ",short_s)), fontface="bold", size=2) +
    geom_tippoint(aes(fill=.data[[grp_col]]), pch=21, alpha=0.9, size=2) + 
    # geom_point(aes(x=branch, y=y+offset_text, fill1=gCF),
    #            data=. %>% filter(!is.na(gCF)),
    #            color="black", pch=21, size=5) + 
    # geom_text(aes(x=branch, y=y+offset_text, label=round(gCF, 1)), color="white", size=1.5) + 
    # geom_point(aes(x=branch, y=y-offset_text, fill2=sCF), 
    #            data=. %>% filter(!is.na(sCF)),
    #            color="black", pch=21, size=5) + 
    # geom_text(aes(x=branch, y=y-offset_text, label=round(sCF, 1)), color="white", size=1.5) + 
    # geom_point(aes(x=branch,y=y, fill3=support), 
    #            data=. %>% filter(!is.na(support), support<=100),
    #            color="black", pch=21, size=5) +
    # geom_text(aes(x=branch, label=round(support, 1)), color="white", size=1.5) + 
    geom_text2(aes(x=branch, label=paste0(support, "/", gCF, "/", sCF), color=mean_sup),
              data = .  %>% filter(!is.na(support), !is.na(gCF), !is.na(sCF)),
              size=2, vjust=-.3, show.legend = F) + 
    geom_treescale() +
    scale_fill_manual(values = col_palette) + 
    scale_color_distiller(palette="RdPu", direction = 1, limits = c(0,100)) + 
    coord_cartesian(clip = 'off') + 
    # scale_x_continuous(expand = c(.2,.2)) +
    guides(fill=guide_legend(title=NULL, nrow=ceiling(length(col_palette)/3))) +
    theme(legend.position="bottom") + # plot.margin=margin(1, 120, 1, 1))
    xlim(c(0, right_x))

ggsave(paste0(prefix, "_consensus.pdf"), cons_plot, units="in", width=7, height=height_plot)

pdf(file = paste0(prefix, "_densitree.pdf"), width=10+right_x, height=height_plot)
par(mar = c(0, 0, 0, 5))
phangorn::densiTree(trees, scaleX = T, scale.bar = F, width=2, consensus = consensus_nms, font = 2, cex = 1)
dev.off()

# Finally save species tree
write.tree(consensus, opt$output)
