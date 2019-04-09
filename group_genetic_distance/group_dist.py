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
To derive within- and between-groups genetic distance
based on pairwise genetic distance matrix, which can be generated with
R ape package `dist.dna(a, as.matrix = TRUE)`
    '''

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-msa_phylip', metavar='<file>', required=False,
        help='the msa file in phylip format')

    parser.add_argument('-pairwise_dist', metavar='<file>', required=False,
        help='file of pairwise genetic distance matrix')

    parser.add_argument('-group_definition', metavar='<file>', required=True,
        help='File format: seqid groupName')

    parser.add_argument('-delimiter', metavar='<str>', default=r'\s+',
        help='the delimiter between `seqid groupName` [%(default)s]')

    parser.add_argument('-b_o', metavar='<file>', type=argparse.FileType('w'),
        default=sys.stdout, help='between-groups distance output')

    parser.add_argument('-i_o', metavar='<file>', type=argparse.FileType('w'),
        default=sys.stdout, help='within-group distance output')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if not args.pairwise_dist and not args.msa_phylip:
        sys.exit('You must specify either `-msa_phylip` or `-pairwise_dist`!')

    return args


def write_R_tmp_script():
    script = '''
require(ape)

myArgs <- commandArgs(trailingOnly = TRUE)
a = read.dna(myArgs[1])
M <- dist.dna(a, as.matrix = TRUE)

write.csv(M, file=myArgs[2])
    '''

    tmp_file = 'tmp_get_dist_pairwise_matrix.R'
    with open(tmp_file, 'w') as fhout:
        fhout.write(script)

    return tmp_file


def get_dist_pairwise_matrix(msa_phylip=None, R_script=None):
    matrix_file = os.path.basename(msa_phylip) + '.csv'
    cmd = 'Rscript ' + R_script + ' ' + msa_phylip + ' ' +  matrix_file
    subprocess.check_output(cmd, shell=True)

    print('matrix file: ', matrix_file)

    return matrix_file


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
        return 'NA.', 'NA.', 'NA.'
    return min(vals), max(vals), statistics.mean(vals)


def get_members_of_paired_group(mdict=None, group_seqid=None, group1=None, group2=None):
    items = []
    for key1 in mdict:
        for key2, val in mdict[key1].items():
            if key1 in group_seqid[group1] and key2 in group_seqid[group2]:
                items.append((key1, key2, val))
                continue
            if key2 in group_seqid[group1] and key1 in group_seqid[group2]:
                items.append((key1, key2, val))
                continue
    return items


def dist_within_group(mdict=None, group_seqid=None, out_handle=None):
    print('# within-group distances:', file=out_handle)
    print('# group\tMinimum\tMaximum\tAverage', file=out_handle)
    for group in sorted(group_seqid.keys()):
        items = get_members_of_specific_group(mdict=mdict, group_seqid=group_seqid, group=group)
        line = [str(i) for i in stat_dist(items)]
        print(group, '\t'.join(line), sep='\t', file=out_handle)


def dist_between_groups(mdict=None, group_seqid=None, out_handle=None):
    groups = sorted(group_seqid.keys())
    already_computed_pairs = {}
    print('# between-groups distances:', file=out_handle)
    print('# group1\tgroup2\tMinimum\tMaximum\tAverage', file=out_handle)
    for group1 in groups:
        for group2 in groups:
            if group1 == group2:
                continue
            forward = '{0} : {1} '.format(group1, group2)
            backward = '{0} : {1} '.format(group2, group1)
            if forward in already_computed_pairs or backward in already_computed_pairs:
                continue
            already_computed_pairs[forward] = 1
            already_computed_pairs[backward] = 1
            items = get_members_of_paired_group(mdict=mdict, group_seqid=group_seqid, group1=group1, group2=group2)
            line = [str(i) for i in stat_dist(items)]
            print(group1, group2, '\t'.join(line), sep='\t', file=out_handle)


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



def main():

    f_dir = os.getcwd()
    R_script = write_R_tmp_script()

    args = get_para()

    if args.msa_phylip:
        args.pairwise_dist = get_dist_pairwise_matrix(
            msa_phylip=args.msa_phylip,
            R_script=R_script)

    triu_dict, tril_dict = csv2dict(args.pairwise_dist, header=0, all_key_to_all=True)

    seqid_group, group_seqid = get_seqid_group(infile=args.group_definition, delimiter=args.delimiter)


    filtered_group_seqid = filter_not_existing_groups(triu_dict=triu_dict, seqid_group=seqid_group, group_seqid=group_seqid)

    dist_between_groups(
        mdict=triu_dict,
        group_seqid=filtered_group_seqid,
        out_handle=args.b_o)

    dist_within_group(
        mdict=triu_dict,
        group_seqid=filtered_group_seqid,
        out_handle=args.i_o)


if __name__ == '__main__':
    main()



