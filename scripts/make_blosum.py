#! /usr/bin/env python

import sys
from collections import OrderedDict


def render_golang(matrix):
    result = []
    result.append('{')
    for key1, row in matrix.items():
        if key1 in ('B', 'Z', 'X'):
            continue
        result.append('\t{}: {{'.format(key1))
        for key2, score in row.items():
            if key2 in ('B', 'Z', 'X'):
                continue
            result.append('\t\t{}: {},'.format(key2, score))
        result.append('\t},')
    result.append('}')
    return '\n'.join(result)


def text2matrix(fname):
    matrix = OrderedDict()
    with open(fname) as fp:
        [*keys] = fp.readline().strip()
        for key1, line in zip(keys, fp):
            vals = line.split()
            for key2, val in zip(keys, vals):
                matrix.setdefault(key1, OrderedDict())[key2] = val
                matrix.setdefault(key2, OrderedDict())[key1] = val
    return matrix


def main():
    fname = sys.argv[1]
    matrix = text2matrix(fname)
    print(render_golang(matrix))


if __name__ == '__main__':
    main()
