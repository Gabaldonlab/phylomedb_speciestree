library(ComplexHeatmap)

# get bl distance matrix
# cat output/data/single_copy/best_trees_0003_sc_labs.nwk | nw_distance -n -m m - > test.txt
# get topological distance matrix 
# cat output/data/single_copy/best_trees_0003_sc_labs.nwk | nw_topology - | nw_distance -n -m m - > test_internode.txt

lines <- readLines("test.txt")

idxs <- which(grepl("^\t", lines))

sps <- unique(unlist(strsplit(paste(lines[idxs], collapse = "\t"), split = "\t")))
sps <- sps[-which(sps=="")]

# sometimes names are only numbers and this fucks up the whole matching
num_sp <- which(grepl('[0-9]', sps))
weird_ones <- sps[num_sp]
sps[num_sp] <- paste0("X", sps[num_sp])


occurence_matrix <- matrix(0, length(sps), length(sps), dimnames = list(sps,sps))
value_matrix <- matrix(0, length(sps), length(sps), dimnames = list(sps,sps))

start <- 1

while (start < length(idxs)){
  
  first_line <- idxs[start]
  last_line <- idxs[start+1]-1
  x <- as.matrix(read.delim(con <- textConnection(lines[first_line:last_line]), row.names = 1))
  close(con)
  n_sp <- ncol(x)
  idxs_weird <- which(rownames(x) %in% weird_ones)
  
  if (length(idxs_weird)>0) {
    
      rownames(x)[idxs_weird] <- paste0("X", rownames(x)[idxs_weird])
      # colnames(x)[idxs_weird] <- paste0(colnames(x)[idxs_weird], "_sp")
  }
  
  start <- start + 1
  
  # add to occurrence
  x2 <- x
  x2[x2>0] <- 1
  diag(x2) <- 1

  mcol <- match(colnames(x),colnames(occurence_matrix))
  mrow <- match(rownames(x),rownames(occurence_matrix))
  occurence_matrix[mrow,mcol] <- occurence_matrix[mrow,mcol]+x2
  
  # value matrix
  # normalize as some genes may have very different total lengths
  x <- x/n_sp
  value_matrix[mrow,mcol] <- value_matrix[mrow,mcol]+x
}


test <- value_matrix/occurence_matrix

Heatmap(test)

