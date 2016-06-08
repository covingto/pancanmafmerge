

import os, os.path, sys
import json, csv

def ndepth(record):
    try:
        if float(record['n_depth']) < 8:
            return 'ndp'
    except:
        pass
    return None

def varkey(record):
    return '|'.join([record[k] for k in ('Tumor_Sample_Barcode', 'Chromosome', 'Start_Position', 
                                         'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')])

def main(args):
    reader = csv.DictReader(args.INPUT, delimiter = '\t')
    for record in reader:
        _ndepth = ndepth(record)
        if _ndepth:
            args.OUTPUT.write('{0}\t{1}\n'.format(varkey(record), _ndepth))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('INPUT', type = argparse.FileType('r'), help = 'input file or stream')
    parser.add_argument('OUTPUT', type = argparse.FileType('w'), help = 'output file or stream')

    args = parser.parse_args()

    main(args)
