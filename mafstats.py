

import os, os.path, sys
import csv, json
from collections import defaultdict
import pprint

def main(args):
    reader = csv.DictReader(args.INPUT, delimiter = '\t')
    struct = defaultdict(dict)
    ctags = defaultdict(int)
    fstruct = defaultdict(int)
    for i, record in enumerate(reader):
        centers = sorted(list(set(record['CENTERS'].split('|'))))
        filters = record['FILTER'].split(',')
        nctag = '|'.join(centers)
        ctags[nctag] += 1
        if 'PASS' not in filters:
            fstruct[nctag] += 1
        for c in centers:
            for f in filters:
                if f not in struct[c]:
                    struct[c][f] = 0
                struct[c][f] += 1
            if 'total' not in struct[c]:
                struct[c]['total'] = 0
            struct[c]['total'] += 1
        if i % 100000 == 0:
            print "Processed %s records" % str(i)
    json.dump({'tags': ctags, 'struct': struct, 'ftags': fstruct}, args.OUTPUT)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('INPUT', type = argparse.FileType('r'), help = 'input file or stream')
    parser.add_argument('OUTPUT', type = argparse.FileType('w'), help = 'output file or stream')
    args = parser.parse_args()

    main(args)


