
import os, os.path, sys
import logging
import hgsc_vcf
import tempfile, shutil

logger = logging.getLogger()                                                               
logger.setLevel(logging.DEBUG)                                                             
ch = logging.StreamHandler()                                                               
ch.setLevel(logging.DEBUG)                                                                 
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')      
ch.setFormatter(formatter)                                                                 
logger.addHandler(ch)                                                                      
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

class FileSplitter(object):
    def __init__(self, reader, seqdict, nrecs = 10000):
        self.reader = reader
        self.seqdict = seqdict
        self.nrecs = nrecs
        self._fileindex = 0
        self.tmpdir = tempfile.mkdtemp()
    
    def _initialize(self, chrom):
        self._fileindex += 1
        writer = hgsc_vcf.Writer(open(os.path.join(self.tmpdir, '%s.%s.split.vcf' % (chrom, str(self._fileindex))), 'w'), self.reader.header)
        writer.write_header()
        return writer

    def split(self):
        # start to split out the input file
        # output pattern is chr.splitnumber.split.vcf
        # all output files have the input header

        # initialize all 
        outfiles = {c:[] for c in self.seqdict.contigs()}
        outcontainer = {c:self._initialize(c) for c in self.seqdict.contigs()}
        outlines = {c:[] for c in self.seqdict.contigs()}
        for record in self.reader:
            c = record['CHROM']
            if c not in outcontainer:
                continue
            outlines[c].append(record)
                
            if sum([len(v) for v in outlines.values()]) > self.nrecs:
                cc, recs = sorted(outlines.items(), key = lambda x: len(x[1]))[-1]
                logger.info("Max recs reached, dumping %s", cc)
                outfiles[cc].append(outcontainer[cc].fobj.name)
                for rec in sorted(recs, key = lambda x: int(x['POS'])):
                    try:
                        outcontainer[cc].write_record(rec)
                    except:
                        logger.error("Failed to write line %s", rec)
                outcontainer[cc].fobj.close()
                outcontainer[cc] = self._initialize(cc)
                outlines[cc] = []
        for c, recs in outlines.items():
            writer = outcontainer[c]
            outfiles[c].append(writer.fobj.name)
            for rec in sorted(recs, key = lambda x: int(x['POS'])):
                try:
                    writer.write_record(rec)
                except:
                    logger.error("Failed to write line %s", rec)
            writer.fobj.close()
        return outfiles

class FileMerger(object):
    def __init__(self, writer, splitfiles, seqdict):
        self.writer = writer
        self.splitfiles = splitfiles
        self.seqdict = seqdict

    def merge(self):
        for c in self.seqdict.contigs():
            logger.info("Merging %s", c)
            mergefiles = self.splitfiles[c]
            self._merge_contig(mergefiles)
            # merge is complete, remove the files
            for f in mergefiles:
                os.remove(f)

    def position_compare(self, x, y):
        return int(x.peek()['POS']) - int(y.peek()['POS'])

    def _merge_contig(self, mergefiles):
        holder = [hgsc_vcf.Reader(open(f, 'r')) for f in mergefiles]
        for r in holder:
            if r.peek() is None:
                r.take() # will iterate if it can
        while len(holder) > 1:
            logger.info("Sorting")
            holder = sorted([r for r in holder if r.peek() is not None], cmp=self.position_compare)
            if len(holder) < 2:
                break
            h0 = holder[0]
            h1 = holder[1]
            while h0.peek() is not None and self.position_compare(h0, h1) < 1:
                self.writer.write_record(h0.take())
        last_r = holder[0]
        while last_r.peek() is not None:
            self.writer.write_record(last_r.take())

def main(args):
    # get the seqdict
    seqdict = SeqDict(args.seqdict)

    # split the file
    splitter = FileSplitter(hgsc_vcf.Reader(open(args.input, 'r')), seqdict)
    splitfiles = splitter.split()
    splitter.reader.fobj.close()

    # write the new files
    merger = FileMerger(hgsc_vcf.Writer(open(args.output, 'w'), splitter.reader.header), splitfiles, seqdict)
    merger.writer.write_header()
    merger.merge()
    merger.writer.fobj.close()
    
    shutil.rmtree(splitter.tmpdir)

    logger.info("Done")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument('seqdict', type = str, help = 'sequence dict file')
    parser.add_argument('input', type = str, help = 'input file')
    parser.add_argument('output', type = str, help = 'output file')

    args = parser.parse_args()

    main(args)
