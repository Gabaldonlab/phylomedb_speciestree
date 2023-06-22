#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from itertools import combinations_with_replacement
import argparse
import os
from ete3 import Tree


def parse_args():
    parser=argparse.ArgumentParser(description="Get distance matrix of phylome homology")
    parser.add_argument('-i', '--input', dest='input', required=True, type=str,
                    help='Input best trees raw file')
    parser.add_argument('-o', '--output', dest='output', required=True, type=str,
                    help='output file')
    parser.add_argument('-m', '--maxtrees', dest='maxtrees', required=False, type=int,
                    help='stop after these trees, useful for testing')
    args=parser.parse_args()
    return args


def pairwise(X, operation):        
    # Initialise precomputed matrix
    precomputed = np.zeros((X.shape[0], X.shape[0]), dtype='int')
    # Initialise iterator over objects in X
    iterator    = combinations_with_replacement(range(X.shape[0]), 2)
    # Perform the operation on each pair
    for i, j in iterator:
        precomputed[i, j] = operation(X[i], X[j])           
    # Make symmetric and return
    return precomputed + precomputed.T - np.diag(np.diag(precomputed))

# define function to compare the files
def overlap(x, y):
    # den=min(len(set(x)), len(set(y)))
    return len(set(x) & set(y))#/den


def read_best_trees(filename):
    tree_dict = {}
    with open(filename) as trees:
        for line in trees:
            line = line.rstrip().split()
            tree_dict[line[0]] = ','.join([leaf.name for leaf in Tree(line[3])])
    return tree_dict



def main():
    inputs=parse_args()
    tree_dict = read_best_trees(inputs.input)
    # output = os.path.splitext(inputs.input)[0]+'_distance.tsv'
    
    # first, read the file
    # df = pd.read_table(inputs.input, names=["seed","hits"],  skip_blank_lines=True)
    df = pd.DataFrame.from_dict(tree_dict, orient='index').reset_index()
    df.columns = ['seed', 'hits']
    # some rows may be empty!
    # print(df. shape[0])
    df.dropna(how="all", inplace=True)
    # print(df. shape[0])
    if inputs.maxtrees is not None:
        df = df.head(inputs.maxtrees)
    sets = [set(el.split(",")) for el in df["hits"]]

    X = np.array(sets)

    dist_X = pairwise(X, overlap)

    df_out = pd.DataFrame(dist_X, columns=df["seed"]).set_index(df["seed"])

    df_out.to_csv(inputs.output, sep="\t")



if __name__ == '__main__':
    main()
