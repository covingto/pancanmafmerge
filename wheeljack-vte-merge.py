

import os, os.path, sys
import bz2
import csv
import re
import glob
from collections import defaultdict
import itertools

import logging

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)



re_site = re.compile(r'\[([0-9]+)~([ATCG]+)>([ATCG]+)\]')

def vtetomafkeys(record):
    oc = record['CALLERS']
    subject = record['SUBJECT']
    chrom = record['CHROM']
    for start, stop, ref, alt in octomafpos(oc):
        yield '|'.join([subject, chrom, str(start), str(stop), ref, alt])

def maftomafkey(record):
    return '|'.join([record['Tumor_Sample_Barcode'][:12], 
                     record['Chromosome'],
                     record['Start_Position'],
                     record['End_Position'],
                     record['Reference_Allele'],
                     record['Tumor_Seq_Allele2']
                     ])

def octomafpos(oc):
    for m in re.finditer(re_site, oc):
        # now we have a match, convert to a maf position
        pos = int(m.group(1))
        ref = m.group(2)
        alt = m.group(3)
        yield vcfpos2maf(pos, ref, alt)

def vcfpos2maf(pos, ref, alt):
    # undo the vcf anchoring
    while ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    if len(ref) == len(alt):
        start, stop = pos, pos + len(ref) - 1
    elif len(ref) < len(alt):
        # this is an insertion
        start, stop = (pos if ref else pos - 1, pos + len(ref) - 1 if ref else pos)
        ref = ref if ref else '-'
    elif len(ref) > len(alt):
        # this is a deletion
        start, stop = (pos, pos + len(ref) - 1)
        alt = alt if alt else '-'
    else:
        raise Exception("Impossible logical block")
    return start, stop, ref, alt


VTEVALIDATIONHEADER = [
    'FC_POS',
    'FC_REF',
    'FC_ALT',
    'TUMORFC_REFCOUNT',
    'TUMORFC_ALTCOUNT',
    'TUMORFC_DP',
    'NORMALFC_REFCOUNT',
    'NORMALFC_ALTCOUNT',
    'NORMALFC_DP',
    'TUMORVALFC_REFCOUNT',
    'TUMORVALFC_ALTCOUNT',
    'TUMORVALFC_DP',
    'NORMALVALFC_REFCOUNT',
    'NORMALVALFC_ALTCOUNT',
    'NORMALVALFC_DP'
    ]

HEADERPREFIXMAP = {
    'TUMOR': 'TUMOR',
    'NORMAL': 'NORMAL',
    'TUMOR_VALIDATION': 'TUMORVAL',
    'TUMOR_VALIDATION_0': 'TUMORVAL',
    'NORMAL_VALIDATION': 'NORMALVAL',
    'NORMAL_VALIDATION_0': 'NORMALVAL'
    }

def batchvte(reader):
    firstrecord = reader.next()
    chrom = firstrecord['CHROM']
    pos = firstrecord['POS']
    ref = firstrecord['REF']
    alt = firstrecord['ALT']
    batch = [firstrecord]
    for record in reader:
        _chrom = record['CHROM']
        _pos = record['POS']
        _ref = record['REF']
        _alt = record['ALT']
        # do we need to emit the batch?
        if chrom != _chrom or pos != _pos or ref != _ref or alt != _alt:
            chrom = _chrom
            pos = _pos
            ref = _ref
            alt = _alt
            yield batch
            batch = []
        batch.append(record)
    yield batch # emit the last batch

def constant_factory(value):
    return itertools.repeat(value).next

def mainmap(args):
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '[map] %(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    logger.addHandler(ch)

    if not args.vte:
        vte = []
    else:
        vte = args.vte
    if args.vteglob:
        vte += glob.glob(args.vteglob)
    # map the vte files
    #
    # vte files are emitted as an annotation 
    # and key
    
    for vtefile in vte:
        with bz2.BZ2File(vtefile, 'r') as fi:
            reader = csv.DictReader(fi, delimiter = '\t')
            for vtebatch in batchvte(reader):
                container = defaultdict(constant_factory('.'))
                # load the container
                # set the prefixs to 0 for the count data
                for h in VTEVALIDATIONHEADER:
                    container[h] = 0
                for record in vtebatch:
                    try:
                        if 'TUMOR' in record['SAMPLE']:
                            prefix = 'TUMOR'
                        else:
                            prefix = 'NORMAL'
                        if 'VAL' in record['SAMPLE']:
                            prefix += 'VAL'
                        container[prefix + 'FC_REFCOUNT'] += int(record['REFCOUNT'])
                        container[prefix + 'FC_ALTCOUNT'] += int(record['ALLELECOUNT'])
                        container[prefix + 'FC_DP'] += int(record['DP'])
                    except:
                        logger.error(container)
                        logger.error(prefix)
                        raise
                fr = vtebatch[0]
                container['FC_POS'] = fr['POS']
                container['FC_REF'] = fr['REF']
                container['FC_ALT'] = fr['ALT']
                for mafkey in vtetomafkeys(fr):
                    args.output.write('\t'.join(
                        [mafkey, 'vter'] +
                        [str(container[k]) for k in VTEVALIDATIONHEADER]
                        ) + '\n')
        print "processed ", vtefile
    # map the maf file
    for maffile in args.maf:
        mafreader = csv.DictReader(maffile, delimiter = '\t')
        for i, record in enumerate(mafreader):
            mafkey = maftomafkey(record)
            args.output.write('\t'.join(
                [mafkey, 'mafr'] + 
                [record[k] for k in mafreader.fieldnames]
                ) + '\n')
            if i % 100000 == 0:
                print "processed", i, "mafrecords"


def mainreduce(args):
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '[reduce] %(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    logger.addHandler(ch)

    # get the MAF header
    mafheader = [c.strip() for c in args.orgmaf.readline().split('\t')]
    outputheader = mafheader + VTEVALIDATIONHEADER
    writer = csv.DictWriter(args.output, fieldnames = outputheader, delimiter = '\t', restval='.')

    writer.writeheader()

    def reducebatch(reader):
        lsplit = [c.strip() for c in reader.next().split('\t')]
        key = lsplit[0]
        batch = [lsplit]
        for row in reader:
            lsplit = [c.strip() for c in row.split('\t')]
            if lsplit[0] != key:
                key = lsplit[0]
                yield batch
                batch = []
            batch.append(lsplit)
        yield batch

    for b in reducebatch(args.input):
        vter = [dict(zip(VTEVALIDATIONHEADER, l[2:])) for l in b if l[1] == 'vter']
        mafr = [dict(zip(mafheader, l[2:])) for l in b if l[1] == 'mafr']
        if not mafr:
            # do something to indicate that we missed a maf record
            logger.warning("No maf records for the following vtes...")
            for v in vter:
                logger.warning(v)
        elif not vter:
            for m in mafr:
                writer.writerow(m)
        else:
            for m in mafr:
                for v in vter:
                    o = m.copy()
                    o.update(v)
                    writer.writerow(o)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parsers = parser.add_subparsers()

    map_parser = parsers.add_parser('map', help = 'run the map function')
    map_parser.add_argument('--maf', type = argparse.FileType('r'), nargs='+', help = 'input maf files')
    map_parser.add_argument('--vte', type = str, nargs='+', help = 'input vte files')
    map_parser.add_argument('--vteglob', type = str, help = "vte file glob")
    map_parser.add_argument('--output', default = sys.stdout, type = argparse.FileType('w'), help = 'output map file')
    map_parser.set_defaults(func=mainmap)

    reduce_parser = parsers.add_parser('reduce', help = 'run the reduce function')
    reduce_parser.add_argument('--orgmaf', type = argparse.FileType('r'), help = 'original maf file, used to get the header')
    reduce_parser.add_argument('--input', type = argparse.FileType('r'), default = sys.stdin, help = 'input map file, after sort')
    reduce_parser.add_argument('--output', type = argparse.FileType('w'), default = sys.stdout, help = 'output maf file with merged headers')
    reduce_parser.set_defaults(func=mainreduce)

    args = parser.parse_args()

    args.func(args)
