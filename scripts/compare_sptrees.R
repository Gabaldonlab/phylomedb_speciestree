
suppressMessages(library(optparse))

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="sptrees separated by comma", dest = "input"),
    make_option(c("-c", "--code"), type="character", default=NULL,
                help="code phylomedb", dest = "code"),
    make_option(c('-s', '--sptree'), type = "character", default = NULL, 
                help = "sptree to root", dest = "sptree"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="output species tree", dest = "output")
    # make_option(c("-a","--approximate"), type = "logical", default=FALSE, action = "store_true",
    #             help="if active, reduce matrix by approximate method", dest = "approx")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# suppressMessages(library(TreeDist))
suppressMessages(library(TreeTools))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ape))
# suppressMessages(library(ggtree))
suppressMessages(library(phangorn))

# opt <- NULL
# opt$code <- "0655"
# opt$input<-"output/sptrees/phylome_0655_sptree_ad.nwk,output/sptrees/phylome_0655_sptree_apro.nwk,output/sptrees/phylome_0655_sptree_asteroid.nwk,output/sptrees/phylome_0655_sptree_astral.nwk,output/sptrees/phylome_0655_sptree_dt.nwk,output/sptrees/phylome_0655_sptree_sprax.nwk,output/sptrees/phylome_0655_sptree_wastral.nwk,output/sptrees/phylome_0655_sptree_fmulrfs.nwk"
# opt$output <- "output/results/phylome_0655/phylome_0655_sptree.nwk"
# opt$sptree <- "output/sptrees/phylome_0655_sptree.nwk"
# code <- opt$code

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = T)

# prefix <- paste0(opt$output,"/result_", code)
prefix <- gsub("_rooted_sptree.nwk", "", opt$output)
print(prefix)

tree_files <- unlist(strsplit(opt$input, ","))

trees <- list()

root_file <- paste0(prefix, "_roots.txt")
cat("method\toutgroup", file=root_file, append=FALSE, sep = "\n")

for (file in tree_files){
    
    method <- gsub(".nwk", "", gsub(".*_", "", file))
    tree <- ape::read.tree(file)

    if (method %in% c("dt", "sprax")) {
        first_clades <- Descendants(tree, tree$edge[1], type = "children")[1:2]
        og_node <- first_clades[which.min(sapply(first_clades, function(x) length(unlist(Descendants(tree, x, type="tips")))))]
        ogs <- tree$tip.label[phangorn::Descendants(tree, og_node, type = "tips")[[1]]]
        
        cat(paste(trimws(method), "\t", ogs), file=root_file, append=TRUE, sep="\n")
        
        final_og <- ifelse(method=="sprax", ogs, NA)
    }
    tree <- unroot(tree)
    tree$edge.length <- NULL
    trees[[method]] <- tree
    # reorder(TreeTools::RootTree(tree, og))
}

class(trees) <- 'multiPhylo'

trees <- reorder(RootTree(trees, final_og))

consensus <- reorder(root(read.tree(opt$sptree), final_og, resolve.root = T))
# consensus <- consensus(trees, p=.5, rooted=TRUE)

rfs <- RF.dist(trees, normalize = T)
rf_matrix <- as.matrix(1 - Matrix::forceSymmetric(as.matrix(rfs), uplo="L"))

png(file = paste0(prefix, "_rf.png"), units="cm", width=15, height=15, res=200)
Heatmap(rf_matrix, 
        "RF_methods",
        col = circlize::colorRamp2(seq(0,1,length.out=9), RColorBrewer::brewer.pal(9, "YlGnBu")),
        width = ncol(rf_matrix)*unit(1, "cm"), 
        height = nrow(rf_matrix)*unit(1, "cm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(rf_matrix[i, j])){
                grid.text(sprintf("%1.2f", rf_matrix[i, j]), x, y, 
                          gp = gpar(fontsize = 8, col=ifelse(rf_matrix[i, j]<.75, "black", "white")))
            } else {
                NULL
            }
        })
dev.off()

if (!is.binary(consensus)) {
    paste("Consensus species tree is not binary. Multifurcations will be solved randomly")
    consensus <- multi2di(consensus)
}

png(file = paste0(prefix, "_densitree.png"), units="in", width=10, height=10, res=200)
phangorn::densiTree(trees, scaleX = T, scale.bar = F, consensus = consensus)
dev.off()

png(file = paste0(prefix, "_consensus.png"), units="in", width=10, height=10, res=200)
# consensus$node.label[which(consensus$node.label == 1)] <- NA
plot(consensus, font = 2, cex = 1)
nodelabels(consensus$node.label, adj = -.05, cex = 0.75, col = "red", frame = "n")
dev.off()

# consensus$node.label <- NULL
write.tree(consensus, opt$output)


# library(ggtree)
# 
# ggtree(consensus, layout = "circular") + 
#     # geom_nodelab(aes(label=gsub("/", "\n", label), color="red")) +
#     geom_nodelab(color="red", hjust = 1) +
#     geom_tiplab2() +
#     theme(legend.position = "none")
#     # geom_nodepoint(aes(fill = as.numeric(unlist(strsplit(label, "/", fixed = T)[1]))))

