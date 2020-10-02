#!/usr/bin/env python3
import re
from mglcmdtools import csv2dict
import statistics
import argparse
import sys
import subprocess
import os


def get_para():
    description = '''
To derive within-groups genetic distance of all genes
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
        default=sys.stdout, help='within-group distance output')

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


def get_members_of_specific_group(mdict=None, group_seqid=None, group=None):
    items = []
    for key1 in mdict:
        for key2, val in mdict[key1].items():
            if key1 in group_seqid[group] and key2 in group_seqid[group]:
                items.append((key1, key2, val))
    return items


def stat_dist(items=None):
    vals = []
    for key1, key2, val in items:
        vals.append(val)

    # one val corresponds to two members.
    if len(vals) <= 0:
        return 'NA.', 'NA.', 'NA.', 'NA.', 'NA.'
    return min(vals), max(vals), statistics.mean(vals), statistics.median(vals), statistics.stdev(vals)


def filter_not_existing_groups(triu_dict=None, seqid_group=None, group_seqid=None):
    # delete some seqids
    filtered_group_seqid = {}
    for group in group_seqid:
        for seqid in group_seqid[group]:
            if seqid in triu_dict:
                if group not in filtered_group_seqid:
                    filtered_group_seqid.setdefault(group, [])
                filtered_group_seqid[group].append(seqid)
            else:
                # discard
                pass

    return filtered_group_seqid


def distStat_of_all_genes_of_same_group(pairwise_dist_list=None, group_seqid=None, out_handle=None):
    group_seq1Seq2Dist = {}
    with open(pairwise_dist_list, 'r') as fh:
        for f in fh:
            f = f.strip()
            triu_dict, tril_dict = csv2dict(f, header=0, all_key_to_all=False)
            filtered_group_seqid = filter_not_existing_groups(triu_dict=triu_dict, group_seqid=group_seqid)

            for group in sorted(filtered_group_seqid.keys()):
                if group not in group_seq1Seq2Dist:
                    group_seq1Seq2Dist.setdefault(group, [])
                items = get_members_of_specific_group(mdict=triu_dict, group_seqid=group_seqid, group=group)
                group_seq1Seq2Dist[group].extend(items)

    print('# within-group distances:', file=out_handle)
    print('# group\tMinimum\tMaximum\tAverage\tMedian\tSample_std', file=out_handle)

    for group in sorted(group_seq1Seq2Dist.keys()):
        items = group_seq1Seq2Dist[group]
        if len(items) < 2:
            print(group, 'NA.', sep='\t', file=out_handle)
            continue
        line = [str(i) for i in stat_dist(items)]
        print(group, '\t'.join(line), sep='\t', file=out_handle)




def main():

    f_dir = os.getcwd()

    args = get_para()

    seqid_group, group_seqid = get_seqid_group(
        infile=args.group_definition,
        delimiter=args.delimiter)

    distStat_of_all_genes_of_same_group(
        pairwise_dist_list=args.pairwise_dist_list,
        group_seqid=group_seqid,
        out_handle=args.i_o)



if __name__ == '__main__':
    main()



