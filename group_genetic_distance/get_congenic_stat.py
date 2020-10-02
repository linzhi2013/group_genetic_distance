#!/usr/bin/env python3
import sys




def main():
    usage = '''
To extract the congenic genetic distance only.

The first two column should be like:
Eothenomys_cachinus     Myodes_rutilus

The genus and species must be seperated by '_'.

lines start with '#' will also be output.


python3 {0} <between-group.dist>
'''.format(sys.argv[0])

    if len(sys.argv) != 2:
        sys.exit(usage)

    d_f = sys.argv[1]
    with open(d_f, 'r') as fh:
        for i in fh:
            i = i.strip()
            if i.startswith('#'):
                print(i)
                continue

            col_1, col_2 = i.split()[0:2]
            g_1 = col_1.split('_')[0]
            g_2 = col_2.split('_')[0]
            if g_1 == g_2 and col_1 != col_2 :
                print(i)


if __name__ == '__main__':
    main()