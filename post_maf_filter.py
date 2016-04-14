

import os, os.path, sys
import gzip, csv
import itertools

def main(args):
    with gzip.open(args.INFILE, 'rb') as fi, gzip.open(args.OUTFILE, 'wb') as fo:
        reader = csv.DictReader(itertools.ifilter(lambda x: x[0]!='#', fi), delimiter = '\t')
        fo.write('# version 2.4\n')
        writer = csv.DictWriter(fo, fieldnames = reader.fieldnames, delimiter = '\t')
        writer.writeheader()
        for i, record in enumerate(reader):
            if i % 1000 == 0:
                print "processed %s lines" % str(i)
            if 'PINDEL' in record['CENTERS']:
                if int(record['t_alt_count']) < 3:
                    continue
            writer.writerow(record)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('INFILE', type = str, help = 'input file')
    parser.add_argument('OUTFILE', type = str, help = 'output file')

    args = parser.parse_args()

    main(args)
