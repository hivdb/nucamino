#! /usr/bin/env python
"""
hfindels.py

Find high frequent insertion/deletion positions.

This script accept an amino acids FASTA file containing MSA
(Multiple Sequence Alignment), it will look for indel positions
appeared in high frequency. The script is useful for discovering
potential positional indel bonus/penalty for NucAmino.

"""
from __future__ import print_function

import sys


def fasta_reader(filename):
    with open(filename) as fp:
        header = None
        seq = []
        for line in fp:
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                if seq:
                    yield header, ''.join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            yield header, ''.join(seq)


def main():
    if len(sys.argv) != 4:
        print('Usage: {} <MIN_NUM_INDELS> <REF_HEADER> <FASTA>'
              .format(sys.argv[0]), file=sys.stderr)
        exit(1)
    _, min_num, ref_header, fasta = sys.argv
    min_num = int(min_num)
    refseq = None
    seqs = list(fasta_reader(fasta))
    for header, seq in seqs:
        if header == ref_header:
            refseq = seq
            break
    else:
        print('Can not locate reference {}'
              .format(ref_header), file=sys.stderr)
        exit(2)
    pos = 1
    seqs = [seq for _, seq in seqs]
    total = len(seqs)
    curins = None
    for refaa, posaas in zip(refseq, zip(*seqs)):
        if refaa == '-':
            if not curins:
                curins = [0, [aa != '-' for aa in posaas]]
            curins[0] += 1
            curins[1] = [
                pre or aa != '-' for pre, aa in zip(curins[1], posaas)]
        else:
            if curins:
                numins = curins[1].count(True)
                if numins >= min_num:
                    print('Insertion at {}: {}/{}'
                          .format(pos - 1, numins, total))
                curins = None
            numdel = posaas.count('-')
            if numdel >= min_num:
                print('Deletion at {}: {}/{}'.format(pos, numdel, total))
            pos += 1


if __name__ == '__main__':
    main()
