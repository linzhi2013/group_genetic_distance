#!/usr/bin/env python3
import re
from mglcmdtools import csv2dict,csv2tupe
import statistics
import argparse
import sys
import subprocess
import os
import numpy as np


def get_para():
    description = '''
To derive between-groups genetic distance of all genes
based on individual gene pairwise genetic distance matrix, which
can be generated with R ape package `dist.dna(a, as.matrix = TRUE)`
    '''

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-pairwise_dist_list', metavar='<file>', required=True,
        help='list of files of pairwise genetic distance matrix')

    parser.add_argument('-group_definition', metavar='<file>', required=True,
        help='File format: seqid groupName')

    parser.add_argument('-delimiter', metavar='<str>', default=r'\s+',
        help='the delimiter between `seqid groupName` [%(default)s]')

    parser.add_argument('-i_o', metavar='<file>', type=argparse.FileType('w'),
        default=sys.stdout, help='between-group distance output')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    return args


def get_seqid_group(infile=None, delimiter=r'\s+'):
    seqid_group = {}
    group_seqid = {}
    with open(infile, 'r') as fh:
        for i in fh:
            i = i.strip()
            seqid, group = re.split(delimiter, i, maxsplit=1)
            seqid_group[seqid] = group
            if group not in group_seqid:
                group_seqid.setdefault(group, [])
            group_seqid[group].append(seqid)

    return seqid_group, group_seqid


def stat_dist(vals=None):
    # one val corresponds to two members.
    if len(vals) <= 0:
        return 'NA.', 'NA.', 'NA.', 'NA.', 'NA.', 'NA.'

    return len(vals), min(vals), max(vals), statistics.mean(vals), statistics.median(vals), statistics.stdev(vals)



def distStat_of_all_genes_of_diff_group(pairwise_dist_list=None, seqid_group=None, out_handle=None):
    group1_group2_seq1Seq2Dist = {}
    with open(pairwise_dist_list, 'r') as fh:
        for f in fh:
            f = f.strip()
            triu_tupe, tril_tupe = csv2tupe(f, header=0)
            for a in triu_tupe:
                k1, k2, val = a
                g1 = seqid_group[k1]
                g2 = seqid_group[k2]
                if str(val) == 'inf':
                    val = np.nan

                if k1 not in group1_group2_seq1Seq2Dist:
                    group1_group2_seq1Seq2Dist.setdefault(g1, {})
                if k2 not in group1_group2_seq1Seq2Dist[g1]:
                    group1_group2_seq1Seq2Dist[g1].setdefault(g2, [])

                group1_group2_seq1Seq2Dist[g1][g2].append(val)

    for g1 in sorted(group1_group2_seq1Seq2Dist.keys()):
        for g2 in sorted(group1_group2_seq1Seq2Dist[g1].keys()):
            vals = group1_group2_seq1Seq2Dist[g1][g2]
            line = [str(i) for i in stat_dist(vals)]
            print(g1, g2, '\t'.join(line), sep='\t', file=out_handle)

def main():

    f_dir = os.getcwd()

    args = get_para()

    seqid_group, group_seqid = get_seqid_group(
        infile=args.group_definition,
        delimiter=args.delimiter)

    distStat_of_all_genes_of_diff_group(
        pairwise_dist_list=args.pairwise_dist_list,
        seqid_group=seqid_group,
        out_handle=args.i_o)



if __name__ == '__main__':
    main()



