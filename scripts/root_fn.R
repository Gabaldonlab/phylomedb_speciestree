library(ape)

roots <- read.delim("output/results/phylome_0301/phylome_0301_roots.txt")
sptree <- "output/sptrees/phylome_0301_sptree.nwk"

final_og <- trimws(roots[roots$method=="sprax ", 2])

t <- read.tree(sptree)
plot(t)

tryCatch({
  root(t, node = getMRCA(t, c("99287", "NITMS")), resolve.root = T)
  },
  error = function(e) {
    root(t, node = getMRCA(t, c("ECOLI", "NITMS")), resolve.root = T)
  }
)

plot(reorder(root(t,  final_og, resolve.root = T)))
plot(reorder(root(t,  node=getMRCA(t, c("ECOLI", "99287")), resolve.root = T)))

test <- read.tree("test.nw")
plot(test, show.node.label = T)
plot(reorder(root(test, node = getMRCA(test, c("Bombina", "Tamias")), resolve.root = T)), show.node.label = T)


root_sptree <- function(tree, og) {
  if (length(og)==1) {
    print("Only one outgroup")
    outree <- root(tree, og, resolve.root=T))
  } else {
    root <- tree$edge[1,1]
    mrca <- getMRCA(tree, og)
    if (root==mrca){
      
    } else {
      if (is.monophyletic(t, final_og)) {
        print("Outgroup is monophyletic!")
        outree <- root(tree, og, resolve.root=T)
      } else {
        print("Outgroup is not monophyletic! Rooting with MRCA")
        outree <- root(tree, node=getMRCA(t, og), resolve.root=T)
      }
    }
  }
  return(reorder(outree))
}

og <- match(final_og, t$tip.label)
n <- length(t$tip.label)
t$tip.label[(1:n)[-og]]

plot(reorder(phytools::reroot(t,  node=getMRCA(t, c("ECOLI", "99287")))))

