

import os, os.path, sys
import csv, gzip
import itertools
from mergesort import mergesort

class SeqDict(object):
    def __init__(self, f):
        self.f = f
        with open(f, 'r') as fi:
            _struct = []
            for line in fi.readlines():
                _struct.append([c.strip() for c in line.split('\t')])
            self._struct = _struct
         
    def contigs(self):
        result = []
        for s in self._struct:
            result.append(s[1].split(':',2)[1])
        return result

def cmprecords(a, b, seqorder):
    achrorder = seqorder.index(a['Chromosome'])
    bchrorder = seqorder.index(b['Chromosome'])
    if achrorder != bchrorder:
        return achrorder - bchrorder
    else:
        return int(a['Start_Position']) - int(b['Start_Position'])


def batch(merges, seqorder):
    def mafkeyfun(myobj):
        class K(object):
            def __init__(self, obj, *args):
                self.obj = obj
            def __cmp__(self, other):
                return cmprecords(self.obj, other.obj, seqorder)

    b = []
    for record in mergesort.mergepeek(merges, mafkeyfun):
        if b and cmprecords(b[0], record, seqorder) != 0:
            yield b
            b = [record] # reset
        else:
            b.append(record)
        # find the lowest record in the merges
    yield b
    yield None # a sentinal value, sends None but does not stop iteration (that happens on the next call)

def record_update(record, r):
    for key in ('Verification_Status', 'Validation_Status', 'Mutation_Status'):
        record_val = record.get(key, None)
        if record_val is None or record_val == '.':
            record[key] = r.get(key, '.')
        elif record_val in ('Unknown', 'Untested', 'Somatic'):
            r_val = r.get(key, None)
            if r_val is not None:
                record[key] = r_val
    # updates in place, this returns None
        
def main(args):
    seqdict = SeqDict(args.SEQDICT)
    contigs = seqdict.contigs()
    primary = csv.DictReader(itertools.ifilter(lambda x: not x.startswith('#'), args.PRIMARY), delimiter = '\t')
    merges = [itertools.ifilter(lambda x: x['Chromosome'] in contigs, csv.DictReader(itertools.ifilter(lambda x: not x.startswith('#'), i), delimiter = '\t')) for i in args.MERGES]
    # assume that merges and primary are sorted

    writer = csv.DictWriter(args.OUTPUT, delimiter = '\t', fieldnames = primary.fieldnames if 'MERGESOURCE' in primary.fieldnames else primary.fieldnames + ['MERGESOURCE'], extrasaction = 'ignore', restval = '.')
    primary = itertools.ifilter(lambda x: x['Chromosome'] in contigs, primary)
    writer.writeheader()
    batchiter = batch(merges, contigs)
    nextbatch = batchiter.next()
    print "Batching %(Chromosome)s|%(Start_Position)s|%(Hugo_Symbol)s" % nextbatch[0]
    record = primary.next()
    while record and nextbatch:
        cmp = cmprecords(record, nextbatch[0], contigs)
        if cmp < 0:
            record['MERGESOURCE'] = 'PRIMARY'
            writer.writerow(record)
            try:
                record = primary.next()
            except StopIteration:
                record = None
        elif cmp > 0:
            # we have passed the batch so write out the remaining items
            for r in nextbatch:
                r['MERGESOURCE'] = 'PASTE'
                writer.writerow(r)
            nextbatch = batchiter.next() # move up the batch records
            if nextbatch:
                print "Batching %(Chromosome)s|%(Start_Position)s|%(Hugo_Symbol)s" % nextbatch[0]
            else:
                print "Done with batches!"
        else:
            # we hit paydirt :)
            newnextbatch = []
            for r in nextbatch:
                for matchstr in ('End_Position', 'Tumor_Sample_Barcode', 'Tumor_Seq_Allele2', 'Reference_Allele'):
                    if r.get(matchstr, None) != record.get(matchstr, None):
                        newnextbatch.append(r)
                        break # breaks the for
                else:
                    # case where all match
                    # this is a match and so we should update the record with the info in r
                    record_update(record, r) # this imprints the value and we can drop the r from the nextbatch entries
            record['MERGESOURCE'] = 'MERGE'
            writer.writerow(record)
            try:
                record = primary.next()
            except StopIteration:
                record = None
            if newnextbatch:
                nextbatch = newnextbatch
            else:
                nextbatch = batchiter.next() # we used up all of the records at this position
            print "Batching %(Chromosome)s|%(Start_Position)s|%(Hugo_Symbol)s" % nextbatch[0]
    if record is not None:
        record['MERGESOURCE'] = 'PRIMARY'
        writer.writerow(record)
        for record in primary: # lets keep reading
            record['MERGESOURCE'] = 'PRIMARY'
            writer.writerow(record)
    if nextbatch is not None:
        for r in nextbatch:
            r['MERGESOURCE'] = 'PASTE'
            writer.writerow(r)
        for nextbatch in batchiter: 
            for r in nextbatch:
                r['MERGESOURCE'] = 'PASTE'
                writer.writerow(r)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('PRIMARY', type = argparse.FileType('r'), help = 'input file or stream')
    parser.add_argument('SEQDICT', type = str, help = 'input sequence dict as generated by Picard')
    parser.add_argument('OUTPUT', type = argparse.FileType('w'), help = 'output file')
    parser.add_argument('MERGES', type = argparse.FileType('r'), nargs='+', help = 'input file or streams to merge')

    args = parser.parse_args()

    main(args)
