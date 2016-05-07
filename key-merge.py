##
# This file merges filter marks into a MAF file.
#
# A mark key is generated for each record in the MAF file.
# The input MAF is then sorted by that mark key.
# Simultaneously input filter files are sorted by their mark keys.
#
# The mark keys are Tumor_Sample_Barcode|Chromosome|Start_Position|End_Position|Reference_Allele|Tumor_Seq_Allele2.
# 
# Sorted files are then streamed into the merger subroutines.  As they are sorted, 
# we do not need to hold all of the data in memory at one time.
#
# Input keys are read in batches such that all records with the same key are grouped.
# The marks for these keys are extracted, split by comma, merged, sorted and joined by comma.
# If the length of the marks for a variant key is zero, PASS is reported.

import os, os.path, sys
import csv
from Queue import Queue
from mergesort import mergesort
import threading


def mafkeyfun(record):
    return '|'.join([record['Tumor_Sample_Barcode'], record['Chromosome'], 
                     record['Start_position'], record['End_position'], 
                     record['Reference_Allele'], record['Tumor_Seq_Allele2']])

def markkeyfun(record):
    return record[0]

def batch(iter, keyfunction):
    lastkey = None
    _batch = []
    for r in iter:
        thiskey = keyfunction(r)
        if lastkey and thiskey != lastkey:
            # print 'yield batch %s' % lastkey
            yield lastkey, _batch
            _batch = []
        lastkey = thiskey
        _batch.append(r)
    yield lastkey, _batch

##
# @param maf - sorted maf reader object
# @param sorted_marks - PeekWrappeer wrapped sorted mark reader objects
def merge(maf, sorted_marks):
    batcher = batch(sorted_marks, markkeyfun)
    mkey, mvals = batcher.next() 
    for r in maf:
        rkey = mafkeyfun(r)
        # print "checking key %s" % rkey
        while mkey is not None and rkey > mkey:
            try:
                mkey, mvals = batcher.next() 
            except StopIteration:
                mkey = mvals = None
        if rkey == mkey:
            # print "keymatch %s" % rkey
            fkeys = []
            for mv in mvals:
                # print mv
                fkeys.extend(mv[1].split(','))
            fkeys = sorted(fkeys)
            r['FILTER'] = ','.join(fkeys)
        else:
            r['FILTER'] = 'PASS'
        yield r 
    
def main(args):
    sort_queue = Queue() # these are files that need to be sorted 
    sorted_marks = []
    
    for markinput in args.MARKFILES:
        sort_queue.put(markinput)
    
    def _t_sort_mark():
        while True:
            thisFile = sort_queue.get()
            try:
                # the PeekWrapper will call next() on init.
                thisSort = mergesort.PeekWrapper(
                                mergesort.mergesort(
                                    csv.reader(open(thisFile, 'r'), delimiter = '\t'),
                                    markkeyfun)
                                )
                sorted_marks.append(thisSort)
                print "Mark sort completed: %s" % thisFile
            finally:
                sort_queue.task_done()
    
    for i in range(4):
        t = threading.Thread(target=_t_sort_mark)
        t.daemon = True 
        t.start()
    
    print "Mark sorting started"
    mafreader = csv.DictReader(open(args.maf, 'r'), delimiter = '\t')
    maf = mergesort.PeekWrapper(mergesort.mergesort(
                mafreader,
                mafkeyfun))
    print "MAF sorting is done"
    sort_queue.join()
    print "Mark sorting complete"
    print "Writing"
    with open(args.output, 'w') as fo:
        fields = mafreader.fieldnames
        if 'FILTER' not in fields:
            fields.append('FILTER')
        writer = csv.DictWriter(fo, fieldnames = fields, delimiter = '\t')
        writer.writeheader()
        for r in merge(maf, mergesort.mergepeek(sorted_marks, markkeyfun)):
            writer.writerow(r)
    print "Done"

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--maf', type = str, help = 'input maf file')
    parser.add_argument('--output', type = str, help = 'output maf file')
    parser.add_argument('MARKFILES', nargs = '+', help = 'input mark files')
    
    args = parser.parse_args()
    
    if not args.maf:
        print "--maf is a required argument"
        sys.exit(2)
    if not args.output:
        print '--output is a required argument'
    
    main(args)
    