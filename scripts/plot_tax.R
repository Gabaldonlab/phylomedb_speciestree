library(tidyverse)
library(ape)
library(ggtree)

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="info file", dest = "input"),
  make_option(c("-o", "--output"), type="character", default=NA,
              help="output filename", metavar="character", dest = "output")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)
names_info <- c("taxid", "mnemo", "date", "longest", "source", "sp_og",
                "k","p","c","o","f","g","s")

info <- read_delim(opt$input, col_names = names_info) %>% 
  mutate(mnemo=gsub("\\..*", "", mnemo),
         s=paste0(s, " (",mnemo,")")) %>% 
  column_to_rownames("mnemo") %>% 
  select(k, p, c, o, f, g, s) %>% 
  mutate_all(as.factor)

frm <- ~k/p/c/o/f/g/s
tree <- ape::as.phylo(frm, data = info)

p <- ggtree(tree, color="grey70") + 
  geom_tiplab(color="grey20") +
  geom_nodelab(color="#009E73", hjust = 0) +
  xlim(c(0,6))

ggsave(opt$output, p, width = 7, height = 10, dpi = 300)