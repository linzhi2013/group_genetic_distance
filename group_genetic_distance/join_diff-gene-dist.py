#!/usr/bin/env python3
import sys
import collections



def read_within_group_dist(f_lst=None):
    '''
    The first column will be the key

    per-file-content-format:
    # group Minimum Maximum Average
    '''

    group_vals = {}
    with open(f_lst, 'r') as fh:
        for f in fh:
            f = f.strip()
            group_vals.setdefault(f, {})
            with open(f, 'r') as fh_dist:
                for i in fh_dist:
                    i = i.strip()
                    if i.startswith('#'):
                        continue
                    line = i.split()
                    group = line[0]
                    val = '\t'.join(line[1:])
                    group_vals[f][group] = val

    return group_vals


def read_between_groups_dist(f_lst=None):
    '''
    The first and two columns will be the key

    per-file-content-format:
    # group1 group2 Minimum Maximum Average
    '''

    group_vals = {}
    with open(f_lst, 'r') as fh:
        for f in fh:
            f = f.strip()
            group_vals.setdefault(f, {})
            with open(f, 'r') as fh_dist:
                for i in fh_dist:
                    i = i.strip()
                    if i.startswith('#'):
                        continue
                    line = i.split()
                    group = '---'.join(line[0:2])
                    val = '\t'.join(line[2:])
                    group_vals[f][group] = val

    return group_vals


def join_dist(group_vals=None, val_width=3):
    # first, get a combined group names
    all_groups = set()
    for f in group_vals:
        for group in group_vals[f]:
            all_groups.add(group)

    all_groups = sorted(all_groups)

    # second, create the joined groups
    joined_group_vals = collections.OrderedDict()
    for group in all_groups:
        joined_group_vals.setdefault(group, [])

    f_order = sorted(group_vals.keys())
    for f in f_order:
        for group in all_groups:
            vals = ['NA.'] * val_width
            vals = '\t'.join(vals)
            if group in group_vals[f]:
                vals = group_vals[f][group]
            joined_group_vals[group].append(vals)


    tmp_g = list(all_groups)[0]
    if '---' in tmp_g:
        title = ['#group1', 'group2']
    else:
        title = ['#group']

    for f in f_order:
        title.append(f)
        vals = [''] * (val_width-1)
        title.append('\t'.join(vals))
    title = '\t'.join(title)
    print(title)
    for group in joined_group_vals:
        line = '\t'.join(joined_group_vals[group])
        if '---' in group:
            g1, g2 = group.split('---')
            print(g1, g2, line, sep='\t')
        else:
            print(group, line, sep='\t')


def main():
    usage = '''
python3 {0} <b|w> <f_lst> <val_width>

b: input is between groups data
w: input is within group data

    '''.format(sys.argv[0])

    if len(sys.argv) != 4:
        sys.exit(usage)

    d_type, f_lst, val_width = sys.argv[1:4]

    val_width = int(val_width)

    if d_type == 'b':
        group_vals = read_between_groups_dist(f_lst=f_lst)
    else:
        group_vals = read_within_group_dist(f_lst=f_lst)

    join_dist(group_vals=group_vals, val_width=3)


if __name__ == '__main__':
    main()