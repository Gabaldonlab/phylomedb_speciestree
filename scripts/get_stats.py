#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True)
# parser.add_argument("-s", "--sptree", type=str, required=True)
# parser.add_argument("-a", "--alndir", type=str, required=True)
parser.add_argument("-n", "--numsp", type=int, required=True)
parser.add_argument("-o", "--output", type=str, default="tree_props.csv")

from ete3 import PhyloTree
import pandas as pd
import numpy as np


def get_species(node):
    '''
    Get species name

    Args:
        node (TreeNode): tree node with sequence name

    Returns:
        string: the tree node species label
    '''

    if '_' in node:
        return node.split("_")[1]
    else:
        return node


def read_treeline(line):
    '''
    Function for reading a phylome tree line

    The function gets the line and splits it using the tab delimiter, then
    stores the tree seed, the inference model used, its likelihood and the
    tree as an ete3 PhyloTree object in a dictionary.

    Args:
        line (str): tree line

    Returns:
        dict: containing the seed, the model, the likelihood and the tree

    Raises:
        KeyError: there is not all the information in the line

    '''

    line = line.split('\t')
    tree_dict = dict()

    try:
        tree_dict['seed'] = line[0]
        tree_dict['model'] = line[1]
        tree_dict['likelihood'] = line[2]
        tree_dict['tree'] = PhyloTree(line[3],
                                           sp_naming_function=get_species)
    except KeyError:
        raise 'The tree line does not contain all the information.'

    return tree_dict


def tree_stats_deep(tree, numsp):
    '''
    Get various tree stats

    From a tree the function retrieves basic numerical information about the
    branches lengths. First, it creates a list where stores all the root to
    tip distances, then calculates the median, the mean, the width and the sum

    Args:
        tree (PhyloTree): ete3 PhyloTree object with a get_species_tag function

    Returns:
        dictionary: dictionary with the tree statistics

    Raises:
        Exception: description
    '''

    # Getting all the leaves
    ndlf = tree.get_leaves()

    # Retrieving the root to tip distances
    distl = list()
    for leaf in ndlf:
        distl.append(tree.get_distance(leaf))

    total_dist = sum([node.dist for node in tree.traverse()])
    internal_dist = sum([node.dist for node in tree.traverse() if not node.is_leaf()])
    # get Mean support values
    support = list()
    for node in tree.traverse():
        if not node.is_leaf():
            support.append(node.support)

    nms = tree.get_leaf_names()

    # patristic distance
    # this is very very slow but could be interesting
    # pat_distance = [tree.get_distance(a, b) for idx, a in enumerate(nms) for b in nms[idx+1:]]
    # mean_pat = np.mean(pat_distance)

    # Get ndups and n speciation trees
    # ntrees, ndups, iterator = tree.get_speciation_trees()

    # Generating the output dictionary
    nodedict = dict()
    nodedict['leafno'] = len(distl)
    nodedict['spno'] = len(tree.get_species())
    nodedict['mean_r2t'] = np.mean(distl)
    # variance of root to tip length
    nodedict['var_r2t'] = np.var(distl)
    # max leaf length
    nodedict['width'] = tree.get_farthest_leaf()[1]
    # total tree length
    nodedict['tree_length'] = total_dist
    nodedict['rate'] = nodedict['tree_length']/nodedict['leafno']
    nodedict['mean_bs'] = np.mean(support)
    if nodedict['tree_length']==0:
        nodedict['treeness'] = "NA"
    else:
        nodedict['treeness'] = internal_dist/nodedict['tree_length']
    nodedict['dup_index'] = nodedict['leafno']/nodedict['spno']

    nodedict['occupancy'] = nodedict['spno']/numsp

    if nodedict['spno']==nodedict['leafno']:
        nodedict['single_copy'] = True
    else:
        nodedict['single_copy'] = False
    # nodedict['ntrees'] = ntrees
    # nodedict['ndups'] = ndups
    # nodedict['mean_pat'] = mean_pat
    return nodedict


args = parser.parse_args()

gene_file = args.input
# algs_dir = args.alndir
# sp_file = args.sptree

if __name__ == "__main__":
    with open(gene_file) as input:
        genes_dict = [read_treeline(line) for line in input]
    print('Done loading genes!')

    # out_tree_path = '%s.trees.tsv' % (args.outfile)
    # out_aln_path = '%s.aln.tsv' % (args.outfile)

    dict_tree = {}
    # aln_results = []

    # optimize alignment code, it is very slow
    # maybe also parallelize this for loop
    for line in genes_dict:
        tree = line["tree"]
        # tree.set_outgroup(fn.root(tree, sp2age))
        dict_tree[line["seed"]] = tree_stats_deep(tree, args.numsp)
        
        # for mode in ["raw", "clean"]:
        #     align = read_aln(line["seed"], algs_dir, mode=mode)
        #     dict_aln = get_aln_stats(align)
        #     dict_aln["gene"] = line["seed"]
        #     dict_aln["mode"] = mode
        #     aln_results.append(dict_aln)
            # dict_tree[line["seed"]].update(dict_aln[line["seed"]])


    out_df_tree = pd.DataFrame.from_dict(dict_tree, orient='index')
    out_df_tree.to_csv(args.output, index_label="gene", sep="\t", float_format='%11.4f')

    # out_df_aln = pd.DataFrame(aln_results)
    # out_df_aln.to_csv(out_aln_path, index=False, sep="\t")
