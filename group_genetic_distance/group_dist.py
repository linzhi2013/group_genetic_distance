#!/usr/bin/env python3
import re
from mglcmdtools import csv2dict
import statistics
import argparse
import sys
import subprocess
import os
import math

def get_para():
    description = '''
To derive within- and between-groups genetic distance
based on pairwise genetic distance matrix ('-pairwise_dist'), which can be generated with
R ape package `dist.dna(a, as.matrix = TRUE) ('-msa_file')`

Copyright Guanliang Meng (linzhi2012'@'gmail'.'com)
    '''

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-msa_file', metavar='<file>', required=False,
        help='input msa file')

    parser.add_argument('-msa_format', default="fasta",
        choices=["interleaved", "sequential", "clustal", "fasta"],
        help='MSA format [%(default)s]')

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
    
    parser.add_argument('-max_dist', metavar='<float>', type=float, 
        default=10, help="distance greater than this value will be excluded. This could be due to a bug in the dist.dna() function. I observed it happens when the two sequences have long gaps(-) on 5' and 3' respectively while the overlapping region is several bp in length! [%(default)s]")

    parser.add_argument('-model', default="K80",
        choices=["raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84",
"BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock"],
        help='Evolutionary model [%(default)s]')

    parser.add_argument('-pairwise_deletion', action='store_true',
        help='''a logical indicating whether to delete the sites with missing data in a pairwise
way. The default is to delete the sites with at least one missing data for all
sequences (ignored if model = "indel" or "indelblock"). [%(default)s]''')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if not args.pairwise_dist and not args.msa_file:
        sys.exit('You must specify either `-msa_file` or `-pairwise_dist`!')

    return args


def write_R_tmp_script(msa_format='fasta', model="K80", pairwise_deletion=False):
    '''
    https://cran.r-project.org/web/packages/ape/ape.pdf

    • raw, N: This is simply the proportion or the number of sites that differ between each pair of
    sequences. This may be useful to draw “saturation plots”. The options variance and gamma
    have no effect, but pairwise.deletion can.

    • TS, TV: These are the numbers of transitions and transversions, respectively.

    • JC69: This model was developed by Jukes and Cantor (1969). It assumes that all substitutions
    (i.e. a change of a base by another one) have the same probability. This probability is the same
    for all sites along the DNA sequence. This last assumption can be relaxed by assuming that
    the substition rate varies among site following a gamma distribution which parameter must be
    given by the user. By default, no gamma correction is applied. Another assumption is that the
    base frequencies are balanced and thus equal to 0.25.

    • K80: The distance derived by Kimura (1980), sometimes referred to as “Kimura’s 2-parameters
    distance”, has the same underlying assumptions than the Jukes–Cantor distance except that
    two kinds of substitutions are considered: transitions (A <-> G, C <-> T), and transversions
    (A <-> C, A <-> T, C <-> G, G <-> T). They are assumed to have different probabilities. A
    transition is the substitution of a purine (C, T) by another one, or the substitution of a pyrimidine (A, G) by another one. A transversion is the substitution of a purine by a pyrimidine,
    or vice-versa. Both transition and transversion rates are the same for all sites along the DNA
    sequence. Jin and Nei (1990) modified the Kimura model to allow for variation among sites
    following a gamma distribution. Like for the Jukes–Cantor model, the gamma parameter must
    be given by the user. By default, no gamma correction is applied.

    '''
    pairwise_deletion_val = 'FALSE'
    if pairwise_deletion:
        pairwise_deletion_val = 'TRUE'

    script = '''
require(ape)

myArgs <- commandArgs(trailingOnly = TRUE)
a = read.dna(myArgs[1], format="{msa_format}")
M <- dist.dna(a, model = "{model}", pairwise.deletion = {pairwise_deletion_val}, as.matrix = TRUE)

write.csv(M, file=myArgs[2])
    '''.format(msa_format=msa_format, model=model, pairwise_deletion_val=pairwise_deletion_val)

    tmp_file = 'tmp_get_dist_pairwise_matrix.R'
    if not os.path.exists(tmp_file):
        with open(tmp_file, 'w') as fhout:
            fhout.write(script)
    
    return tmp_file


def get_dist_pairwise_matrix(msa_file=None, R_script=None):
    matrix_file = os.path.basename(msa_file) + '.csv'
    cmd = 'Rscript ' + R_script + ' ' + msa_file + ' ' +  matrix_file
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


def stat_dist(items=None, max_dist=10):
    raw_vals = []
    for key1, key2, val in items:
        raw_vals.append(val)

    # remove potential inf
    raw2_vals = [x for x in raw_vals if not math.isinf(x)]
    vals = []
    for x in raw2_vals:
        if x > max_dist:
            print("excluding", x) 
        else:
           vals.append(x)
    # [x for x in vals if x < max_dist]
    # one val corresponds to two members.
    if len(vals) <= 0:
        return 'NA.', 'NA.', 'NA.', 'NA.', 'NA.'
    return min(vals), max(vals), statistics.mean(vals), statistics.median(vals), statistics.stdev(vals)


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


def dist_within_group(mdict=None, group_seqid=None, out_handle=None, max_dist=10):
    print('# within-group distances:', file=out_handle)
    print('# group\tMinimum\tMaximum\tAverage\tMedian\tSample_std', file=out_handle)
    for group in sorted(group_seqid.keys()):
        items = get_members_of_specific_group(mdict=mdict, group_seqid=group_seqid, group=group)
        line = [str(i) for i in stat_dist(items, max_dist)]
        print(group, '\t'.join(line), sep='\t', file=out_handle)


def dist_between_groups(mdict=None, group_seqid=None, out_handle=None, max_dist=10):
    groups = sorted(group_seqid.keys())
    already_computed_pairs = {}
    print('# between-groups distances:', file=out_handle)
    print('# group1\tgroup2\tMinimum\tMaximum\tAverage\tMedian\tSample_std', file=out_handle)
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
            line = [str(i) for i in stat_dist(items, max_dist)]
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
    args = get_para()

    f_dir = os.getcwd()
    R_script = write_R_tmp_script(
        msa_format=args.msa_format,
        model = args.model,
        pairwise_deletion=args.pairwise_deletion)


    if args.msa_file:
        args.pairwise_dist = get_dist_pairwise_matrix(
            msa_file=args.msa_file,
            R_script=R_script)

    triu_dict, tril_dict = csv2dict(args.pairwise_dist, header=0, all_key_to_all=True)

    seqid_group, group_seqid = get_seqid_group(infile=args.group_definition, delimiter=args.delimiter)


    filtered_group_seqid = filter_not_existing_groups(triu_dict=triu_dict, seqid_group=seqid_group, group_seqid=group_seqid)

    dist_between_groups(
        mdict=triu_dict,
        group_seqid=filtered_group_seqid,
        out_handle=args.b_o,
        max_dist=args.max_dist)

    dist_within_group(
        mdict=triu_dict,
        group_seqid=filtered_group_seqid,
        out_handle=args.i_o,
        max_dist=args.max_dist)


if __name__ == '__main__':
    main()



