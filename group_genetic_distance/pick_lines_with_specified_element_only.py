#!/usr/bin/env python3
import sys
import os
import re


def check_match(c, elements=None):
    for e in elements:
        if e in c:
            return True

    return False


def main():
    usage = '''
python3 {} <infile> <e1> <e2> [e3] [...]

To output the lines whose first two columns consist of any combinations of user specified elements.

NB: The match is checked within the 'in' function!

'''.format(sys.argv[0])

    if len(sys.argv) < 3:
        sys.exit(usage)

    in_f = sys.argv[1]
    elements = sys.argv[2:]


    with open(in_f, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i:
                continue

            c1, c2 = i.split("\t")[0:2]
            if check_match(c1, elements) and check_match(c2, elements):
                print(i)

if __name__ == '__main__':
    main()