#!/usr/bin/env python
"""
This file is contains code to be used alongside pplacer as described in the upcoming ALCOB conference paper :

Wedell, E., Cai, Y., Warnow, T. (2021). Scalable and Accurate Phylogenetic Placementusing pplacer-XR.

Copyright (c) 2021 pplacer-SCAMPP Developers
Yirong Cai <yirongc2@illinois.edu>
Eleanor Wedell <ewedell2@illinois.edu>
All rights reserved.

Licence: MIT Licence, 
see https://opensource.org/licenses/MIT

pplacer can be found at: 
https://github.com/matsen/pplacer

***must be run from the directory containing pplacer***
"""

import sys
import os
import utils
import shutil
import json
import time
import argparse
import treeswift

def main(args):
    tree_path = args.tree
    output = args.outdir
    outFile = args.output
    aln = args.alignment
    refaln = args.refaln
    n = args.subtreesize
    run = args.tmpfilenbr
    subtree_flag = args.subtreetype
    frag_flag = args.fragmentflag
    q_aln = args.qalignment
    model = args.model
    info = args.info

    # read msa and reference tree
    t0 = time.perf_counter()
    tree = treeswift.read_tree_newick(tree_path)

    leaf_dict = tree.label_to_node(selection=set([n.label for n in tree.traverse_leaves()])) #selection='leaves')

    for leaf in leaf_dict:
        assert leaf != ''
    
    if tree.root.num_children() == 2:
        tree.deroot()

    if q_aln != "":
        ref_dict = utils.read_data(aln)
        q_dict = utils.read_data(q_aln)
        #print("ref_dict", ref_dict)
        #print("q_dict", q_dict)
    else:
        aln_dict = utils.read_data(aln)
        #print("reading in from alignment")
        #print("aln_dict", aln_dict)
        ref_dict, q_dict = utils.seperate(aln_dict, leaf_dict)
        #print("ref_dict", ref_dict)
        #print("q_dict", q_dict)
    
    jplace = dict()
    placements = []

    # create edge tokens to be used with jplace output
    utils.add_edge_nbrs(tree)
    jplace["tree"] = utils.newick_edge_tokens(tree)
    print(f'{time.perf_counter() - t0} seconds loading data')

    files = []
    try:
        os.makedirs(f"tmp{run}", exist_ok=True)
    except OSError as error:
        pass
    try:
        os.makedirs(output, exist_ok=True)
    except OSError as error:
        pass

    # place each query sequence
    for name, seq in q_dict.items():
        tmp_tree = f"tmp{run}/tree_{name}"
        tmp_aln = f"tmp{run}/aln{name}.fa"
        tmp_output = f"tmp{run}/{name}.jplace"

        # finds closest sister taxon and subtree leaves
        if subtree_flag == 'h':
            y = utils.find_closest_hamming(seq, ref_dict, n, frag_flag)
            labels = []
            for taxon in y:
               labels.append(leaf_dict[taxon].get_label())
            print(f'Closest sister taxon found: {y[0]}')
        else:
            y = utils.find_closest_hamming(seq, ref_dict, 1, frag_flag)
            print(f'Closest sister taxon found: {y[0]}')
            print (f'{time.perf_counter() - t0} seconds finding closest leaf')
            node_y = leaf_dict[y[0]]
            if subtree_flag == 'n':
                labels = utils.subtree_nodes(tree, node_y, n)
            else:
                labels = utils.subtree_nodes_with_edge_length(tree, node_y, n)

        # write subtree MSA and aligned query sequence to tmp file
        with open(tmp_aln, "w") as f:
            f.write(">"+name)
            f.write("\n")
            f.write(seq+"\n")
            for label in labels:
                label_list = label.split('%%',1)
                label = label_list[0]
                f.write(">"+label+"\n")
                f.write(ref_dict[label])
                f.write("\n")
        subtree = tree.extract_tree_with(labels)

        if subtree.root.num_children() == 2:
            subtree.deroot()
        utils.remove_edge_nbrs(subtree)

        subtree.write_tree_newick(tmp_tree, hide_rooted_prefix=True)

        print(f'{time.perf_counter() - t0} seconds extracting subtree')
        # run pplacer from directory containing pplacer binaries
        #os.system("./pplacer -m {} -s {} -t {} -o {} {}".format(model, info, tmp_tree, tmp_output, tmp_aln))

        # build ref_pkg with info and tmp_tree and tmp_aln
        ref_pkg = "{}/{}.refpkg".format(output, outFile)
        print(f"taxtastic/taxtastic-env/bin/taxit create -P {ref_pkg} -l {name} --aln-fasta {refaln} --tree-file {tmp_tree} --tree-stats {info}")
        os.system(f"taxtastic/taxtastic-env/bin/taxit create -P {ref_pkg} -l {name} --aln-fasta {refaln} --tree-file {tmp_tree} --tree-stats {info}")
        print(f"pplacer -m {model} -c {ref_pkg} -o {tmp_output} -j 1 {tmp_aln} --timing")
        os.system(f"pplacer -m {model} -c {ref_pkg} -o {tmp_output} -j 1 {tmp_aln} --timing")

        print (f'{time.perf_counter() - t0} seconds running pplacer')

        # load the jplace file and find placements in the original backbone tree
        place_file = open(tmp_output, 'r')
        place_json = json.load(place_file)

        if len(place_json["placements"]) > 0:
            added_tree, edge_dict = utils.read_tree_newick_edge_tokens(place_json["tree"])
            tmp_place = place_json["placements"][0]
            for i in range(len(tmp_place["p"])):
                edge_num = tmp_place["p"][i][1] # edge number in subtree
                edge_distal = tmp_place["p"][i][0] # distal length from parent node

                # find placement edge according to edge number
                right_n = edge_dict[str(edge_num)]
                left_n = right_n.get_parent()

                # obtain a path from leaf left to leaf right containing placement edge through the subtree
                left, path_l = utils.find_closest(left_n, {left_n, right_n})
                right, path_r = utils.find_closest(right_n, {left_n, right_n})

                # obtain the corresponding path in backbone tree
                left = leaf_dict[left.get_label()]
                right = leaf_dict[right.get_label()]
                _, path = utils.find_closest(left, {left}, y=right)

                # find the length of placement along the path from leaf left to leaf right in subtree
                length = sum([x.get_edge_length() for x in path_l])+edge_distal

                # find the target placement edge in backbone tree 
                target_edge = path[-1]
                for j in range(len(path)):
                    length -= path[j].get_edge_length()
                    if length < 0:
                        target_edge = path[j]
                        break

                tmp_place["p"][i][0] = 0

                label = target_edge.get_label()
                [taxon, target_edge_nbr] = label.split('%%',1)
                tmp_place["p"][i][0] = target_edge.get_edge_length()+length
                tmp_place["p"][i][1] = int(target_edge_nbr)

            # append the placement to the output jplace
            placements.append(tmp_place.copy())

        place_file.close()

    # build jplace file
    jplace["placements"] = placements
    jplace["metadata"] = {"invocation": " ".join(sys.argv)}
    jplace["version"] = 3
    jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
            "likelihood", "pendant_length"]

    
    output = open('{}/{}.jplace'.format(output,outFile), 'w') 
    json.dump(jplace, output, sort_keys=True , indent=4)
    output.close()
    print(f'{time.perf_counter() - t0} seconds building jplace')
    shutil.rmtree("tmp{}".format(run))
    
def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--info", type=str,
                        help="Path to model parameters", required=True, default=None)
    parser.add_argument("-t", "--tree", type=str,
                        help="Path to reference tree with estimated branch lengths", required=True, default=None)
    parser.add_argument("-d", "--outdir", type=str,
                        help="Directory path for output", required=True, default=None)
    parser.add_argument("-a", "--alignment", type=str,
                        help="Path for query and reference sequence alignment in fasta format", required=True, default=None)
    parser.add_argument("-r", "--refaln", type=str,
                        help="Path for reference sequence alignment in fasta format", required=True, default=None)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file name", required=False, default="pplacer-SCAMPP")  
    parser.add_argument("-m", "--model", type=str,
                        help="Model used for edge distances",
                        required=False, default="GTR")
    parser.add_argument("-b", "--subtreesize", type=int,
                        help="Integer size of the subtree",
                        required=False, default=2000)
    parser.add_argument("-s", "--subtreetype", type=str,
                        help="d (default) for edge weighted distances, n for node distances, h for hamming distances",
                        required=False, default="d")
    parser.add_argument("-n","--tmpfilenbr", type=int,
                        help="tmp file number",
                        required=False, default=0)
    parser.add_argument("-q", "--qalignment", type=str,
                        help="Path to query sequence alignment in fasta format (ref alignment separate)",
                        required=False, default="")
    parser.add_argument("-f", "--fragmentflag", type=bool,
                        help="boolean, True if queries contain fragments",
                        required=False, default=False)
    parser.add_argument("-v", "--version", action="version", version="2.0.0", help="show the version number and exit")
    return parser.parse_args()

if __name__ == "__main__":
    main(parseArgs())
