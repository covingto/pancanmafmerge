
import hgsc_vcf
import os, os.path, sys
import json, re, logging
import threading
from Queue import Queue

logger = logging.getLogger('wheeljack.filter_alts')
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.DEBUG)

def safe_div(a, b):
    try:
        a = float(a)
        b = float(b)
        return a / b
    except:
        return 0.0

def safe_float(a):
    try:
        return float(a)
    except:
        return 0.0

##
# Do the lambdas select this allele in this sample
#
# The allele in this sample must pass all of the selelction lambdas
def _sample_select_filter(sample, lambdas, i, ref_i):
    for f in lambdas:
        try:
            if not eval(f):
                return False
        except:
            logger.error("Error in processing sample: %s", sample)
            logger.error("Info: i: %s, ref_i: %s", i, ref_i)
            raise
    return True

def _samples_filter(samples, lambdas, i, ref_i):
    for s in samples.values():
        if _sample_select_filter(s, lambdas, i, ref_i):
            return True
    return False

def selection_function(record, lambdas):
    gt = hgsc_vcf.split_gt(record)
    ref_i = hgsc_vcf.ref_index(record, gt)
    samples = record['SAMPLES']
    for i, g in enumerate(gt):
        if i == ref_i:
            continue
        if _samples_filter(samples, lambdas, i, ref_i):
            yield i

def build_hgsc_vcf_select_function(lambdas):
    def _custom_select_function(record):
        for i in selection_function(record, lambdas):
            yield i
    return _custom_select_function

def process_vcf(reader, writer, config, simplify = True):
    _process_select_function = build_hgsc_vcf_select_function(config['samplefilter'])
    writer_q = Queue()
    processor_q = Queue()
    def _writer_function():
        i = 0
        while True:
            out_record = writer_q.get()
            try:
                writer.write_record(out_record)
                i += 1
                if not i % 10000:
                    logger.info("Processed %s lines", i)
            except Exception as inst:
                logger.exception("Exception in record writing: %s", str(inst))
            finally:
                writer_q.task_done()
    
    writer_thread = threading.Thread(target=_writer_function)
    writer_thread.daemon = True
    writer_thread.start()
    
    def _processor_function():
        while True:
            in_record = processor_q.get()
            try:
                for b_record in hgsc_vcf.select_allele(in_record, _process_select_function, simplify):
                    writer_q.put(b_record)
            except Exception as inst:
                logger.exception("Exception in record processing: %s", str(inst))
            finally:
                processor_q.task_done()
    
    for i in xrange(1):
        processor_thread = threading.Thread(target=_processor_function)
        processor_thread.daemon = True
        processor_thread.start()
    
    for record in reader:
        if len(record['ALT']) > 0 and record['ALT'][0] != '.':
            processor_q.put(record)
    
    processor_q.join()
    writer_q.join()
    
def main(args):
    reader = hgsc_vcf.Reader(args.INFILE)
    header = reader.header
    writer = hgsc_vcf.Writer(args.OUTFILE, header)
    writer.header.add_header('##COMMAND=<ID=filter-alts,ARGS="%s">' % re.escape(' '.join(sys.argv)))
    writer.write_header()
    config = json.load(args.CONFIG)
    process_vcf(reader, writer, config)
    logger.info("Done")


if __name__ == '__main__':
    mainlogger = logging.getLogger()
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    mainlogger.addHandler(ch)
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('CONFIG', type = argparse.FileType('r'), help = 'config file')
    parser.add_argument('INFILE', type = argparse.FileType('r'), help = 'input VCF file')
    parser.add_argument('OUTFILE', type = argparse.FileType('w'), help = 'output VCF file')

    args = parser.parse_args()


    main(args)
