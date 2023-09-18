library(tidyverse)
library(ggtree)
library(patchwork)
theme_set(theme_bw())

part_files <- list.files("output/input/pgroot/", full.names = T, pattern = "*part_pgroot.csv")
ad_files <- list.files("output/input/pgroot/", full.names = T, pattern = "*ad_pgroot.csv")

for (idx in c(1:length(part_files))) {
  code <- gsub("_part_pgroot.csv", "", gsub(".*phylome_", "", part_files[idx]))
  ad_file <- gsub("part_pgroot.csv", "ad_pgroot.csv", part_files[idx])

  if (!(file.exists(ad_file))) {
    print(code)
    next
  }
  part <- read_delim(part_files[idx], show_col_types = F)
  ad <- read_delim(ad_file, show_col_types = F)
  
  tree <- ape::read.tree(paste0("output/results/phylome_", code, "/phylome_", code, "_rooted_sptree.nwk"))

  mat <- part %>% 
    pivot_longer(!Partition_id) %>% 
    pivot_wider(names_from = Partition_id) %>% 
    mutate_all(as_factor) %>% 
    column_to_rownames("name")
  p <- ggtree(tree) + 
    geom_tiplab(size=2, align=TRUE, linesize=.5)
  
  p <- gheatmap(p, mat, width=1, offset = .2,
           colnames=FALSE, legend_title="Rooting partition") +
    scale_fill_manual(values = c("black", "grey90"))

  p2 <- mutate(ad, isfirst=Partition_id==1) %>% 
  ggplot(aes(AD, group=Partition_id, color=isfirst)) + 
    facet_wrap(~GeneFam) +
    stat_ecdf()

  p <- p / p2
  ggsave(paste0("test/",code, ".png"), p, dpi=300, width = 8, height = 6)
}
