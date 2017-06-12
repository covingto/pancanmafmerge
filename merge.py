

##
# "wrapper" utility for VCF merge and MAF generation

import os, os.path, sys
import tempfile, re
import hgsc_vcf
import subprocess, traceback, shutil

import smtplib
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart

def notify(message, to = ['kylecovington1@gmail.com'], subject="Error in merge", attachments=None):
    me = 'kylecovington1@gmail.com'
    msg = MIMEMultipart(
        From=me,
        To=', '.join(to),
        Subject=subject
    )
    msg.attach(MIMEText('''
Hello,

%s

Thanks, have a nice day!
''' % message))
    msg['Subject'] = subject
    msg['From'] = me
    msg['To'] = ', '.join(to)
    s = smtplib.SMTP('localhost')
    s.sendmail(me, to, msg.as_string())
    s.close()

import logging



logger = logging.getLogger('merge')
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler(stream=sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

logger.addHandler(ch)

PACKAGEDIR = os.path.dirname(os.path.abspath(__file__))

def filter(fpath, caller, tmpdir):
    outputfpath = os.path.join(tmpdir, os.path.splitext(os.path.basename(fpath))[0] + '.filtered.vcf')
    tmpfile = os.path.join(tmpdir, 'tmpfile.vcf')
    if os.path.isfile(outputfpath):
        logger.info("Skipping filtering because %s exists", outputfpath)
        return outputfpath

    if caller.lower() == 'muse':
        logger.info("Applying MuSE filter to %s", fpath)

        subprocess.check_call('python %(PACKAGE)s/filter_muse.py --level 5 %(input)s %(output)s' % {
            'PACKAGE': PACKAGEDIR,
            'input': fpath,
            'output': tmpfile},
            shell = True)

        shutil.move(tmpfile, outputfpath)
        return outputfpath

    elif caller.lower() == 'radia':
        logger.info("Applying RADIA filter to %s", fpath)

        subprocess.check_call('python %(PACKAGE)s/filter_radia.py %(input)s %(output)s' % {
            'PACKAGE': PACKAGEDIR,
            'input': fpath,
            'output': tmpfile},
            shell = True)

        shutil.move(tmpfile, outputfpath)
        return outputfpath

    elif caller.lower() == "somaticsniper":
        nid, tid, nbar, tbar = getTNids(fpath)
        logger.info("Applying SomaticSniper filter to %s", fpath)

        subprocess.check_call('perl %(PACKAGE)s/vcf2maf/vcf2vcf.pl --add-filter --input-vcf %(input)s --output-vcf %(output)s --vcf-tumor-id %(tid)s --vcf-normal-id %(nid)s' %{
            'PACKAGE': PACKAGEDIR,
            'input': fpath,
            'output': tmpfile,
            'tid': tid,
            'nid': nid},
            shell = True)

        shutil.move(tmpfile, outputfpath)
        return outputfpath

    elif caller.lower() == "varscans":
        nid, tid, nbar, tbar = getTNids(fpath)
        logger.info("Applying VarScan SNP filter to %s", fpath)

        subprocess.check_call('perl %(PACKAGE)s/vcf2maf/vcf2vcf.pl --add-filter --input-vcf %(input)s --output-vcf %(output)s --vcf-tumor-id %(tid)s --vcf-normal-id %(nid)s' %{
            'PACKAGE': PACKAGEDIR,
            'input': fpath,
            'output': tmpfile,
            'tid': tid,
            'nid': nid},
            shell = True)

        shutil.move(tmpfile, outputfpath)
        return outputfpath

    elif caller.lower() == "varscani":
        nid, tid, nbar, tbar = getTNids(fpath)
        logger.info("Applying VarScan INDEL filter to %s", fpath)

        subprocess.check_call('perl %(PACKAGE)s/vcf2maf/vcf2vcf.pl --add-filter --input-vcf %(input)s --output-vcf %(output)s --vcf-tumor-id %(tid)s --vcf-normal-id %(nid)s' %{
            'PACKAGE': PACKAGEDIR,
            'input': fpath,
            'output': tmpfile,
            'tid': tid,
            'nid': nid},
            shell = True)

        shutil.move(tmpfile, outputfpath)
        return outputfpath

    else:
        return fpath

def sort(fpath, tmpdir):
    outputfpath = os.path.join(tmpdir, os.path.splitext(os.path.basename(fpath))[0] + '.sorted.vcf')
    tmpfile = os.path.join(tmpdir, 'tmpfile.vcf')
    if os.path.isfile(outputfpath):
        logger.info("Skipping sort because %s exists", outputfpath)
        return outputfpath

    logger.info("Sorting %s -> %s", fpath, outputfpath)

    subprocess.check_call('python %(PACKAGEDIR)s/vcf-sort.py %(seqdict)s %(input)s %(output)s' % {
        'PACKAGEDIR': PACKAGEDIR,
        'seqdict': '/hgsc_software/cancer-analysis/resources/references/human/hg19/hg19.dict',
        'input': fpath,
        'output': tmpfile},
        shell = True)

    shutil.move(tmpfile, outputfpath)
    return outputfpath

##
# returns a 2 tuple with the names of the genotype fields that correspond to the normal and tumor samples
def getTNids(fpath):
    sample_headers = []
    field_header = None
    with open(fpath, 'r') as fi:
        for line in fi.readlines():
            if line[0] != '#': 
                break
            if '##SAMPLE' == line[:8]:
                sample_headers.append(line)
            if '#CHROM' == line[:6]:
                field_header = [c.strip() for c in line.split('\t')]
    primaryline = [l for l in sample_headers if 'ID=PRIMARY' in l]
    normalline = [l for l in sample_headers if 'ID=NORMAL' in l]
    if primaryline:
        primarysplit = hgsc_vcf.metainfo.ComplexHeaderLine(string = primaryline[0])
    else:
        # check if there is a METASTATIC sample (SKCM for example)
        primaryline = [l for l in sample_headers if 'ID=METASTATIC' in l or 'ID=RECURRANCE' in l]
        if primaryline:
            primarysplit = hgsc_vcf.metainfo.ComplexHeaderLine(string = primaryline[0])
            # we did find a METASTATIC sample, convert the field_header
        else:
            raise ValueError("Could not find ID=PRIMARY or ID=METASTATIC or ID=RECURRANCE in file %s, start again :)" % fpath)
    if normalline:
        normalsplit = hgsc_vcf.metainfo.ComplexHeaderLine(string = normalline[0])
    else:
        raise ValueError("Could not find ID=NORMAL in the file %s, start again :)" % fpath)
    samples = [s for s in field_header[9:] if s in ('NORMAL', 'PRIMARY', 'TUMOR', 'METASTATIC', 'RECURRANCE', normalsplit.fields['SampleTCGABarcode'], primarysplit.fields['SampleTCGABarcode'])]
    if len(samples) < 2:
        raise ValueError("Fewer than 2 samples found: %s" % samples)
    elif len(samples) > 2:
        raise ValueError("More than 2 samples found: %s" % samples)
    # ok so we have the data that we need
    if samples[0] == normalsplit.fields['SampleTCGABarcode']:
        return samples[0], samples[1], normalsplit.fields['SampleTCGABarcode'], primarysplit.fields['SampleTCGABarcode']
    elif samples[0] == 'NORMAL':
        return samples[0], samples[1], normalsplit.fields['SampleTCGABarcode'], primarysplit.fields['SampleTCGABarcode']
    elif samples[1] == normalsplit.fields['SampleTCGABarcode']:
        return samples[1], samples[0], normalsplit.fields['SampleTCGABarcode'], primarysplit.fields['SampleTCGABarcode']
    elif samples[1] == 'NORMAL':
        return samples[1], samples[0], normalsplit.fields['SampleTCGABarcode'], primarysplit.fields['SampleTCGABarcode']
    else:
        raise ValueError("Can't figure out the tumor and normal sample id's in %s" % samples)

def v2v(fpath, tmpdir):
    outputfpath = os.path.join(tmpdir, os.path.splitext(os.path.basename(fpath))[0] + '.v2v.vcf')
    tmpfile = os.path.join(tmpdir, 'tmpfile.vcf')
    if os.path.isfile(outputfpath):
        logger.info("Skipping vcf reduction because %s exists", outputfpath)
        return outputfpath
    logger.info("vcf reduction of %s -> %s", fpath, outputfpath)

    nid, tid, nbar, tbar = getTNids(fpath)

    subprocess.check_call('/hgsc_software/perl/perl-5.16.2/bin/perl %(PACKAGEDIR)s/vcf2maf/vcf2vcf.pl --input-vcf %(input)s --output-vcf %(output)s --vcf-tumor-id %(tid)s --vcf-normal-id %(nid)s' % {
        'PACKAGEDIR': PACKAGEDIR,
        'input': fpath,
        'output': tmpfile,
        'tid': tid,
        'nid': nid},
        shell = True)

    shutil.move(tmpfile, outputfpath)
    return outputfpath

def merge(outfile, mergefiles):
    outputfpath = outfile
    tmpfile = os.path.join(os.path.dirname(outputfpath), 'tmpfile.vcf')
    if os.path.isfile(outputfpath):
        logger.info("Skipping merge because %s exists", outputfpath)
        return outputfpath
    logger.info("Merging %s -> %s", mergefiles, outputfpath)

    subprocess.check_call('python %(PACKAGEDIR)s/vcf-merge.py --keys %(keys)s --output %(output)s %(inputs)s' % {
        'PACKAGEDIR': PACKAGEDIR,
        'keys': ' '.join([c for c, f in mergefiles]),
        'output': tmpfile,
        'inputs': ' '.join([f for c, f in mergefiles])},
        shell = True)

    shutil.move(tmpfile, outputfpath)
    return outputfpath

def annotate(fpath, tmpdir):
    outputfpath = os.path.join(tmpdir, os.path.splitext(os.path.basename(fpath))[0] + '.annotated.vcf')
    tmpfile = os.path.join(tmpdir, 'tmpfile.vcf')

    if os.path.isfile(outputfpath):
        logger.info("Skipping annotation because %s exists", outputfpath)
        return outputfpath
    logger.info("Processing annotation of %s -> %s", fpath, outputfpath)
    vepannotation = os.path.join(tmpdir, 'vepannotate.vcf')
    if os.path.isfile(vepannotation):
        logger.info("Skipping vep annotation because %s exists", vepannotation)
    else:
        logger.info("Processing vep annotation")

        subprocess.check_call('export PERL5LIB=/hgsc_software/cancer-analysis/code/vep-82:/users/covingto/perl5/lib/perl5:$PERL5LIB && export PATH=/hgsc_software/cancer-analysis/code/vep-82/htslib:$PATH && /hgsc_software/perl/perl-5.16.2/bin/perl /hgsc_software/cancer-analysis/code/vep-82/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl --dir /hgsc_software/cancer-analysis/code/vep-82/cache/human/grch37 --format vcf --everything -i %(input)s -o %(output)s --cache --vcf --force_overwrite --check_existing --allow_non_variant --buffer_size 100 --offline --fork 2' % {
            'input': fpath,
            'output': tmpfile},
            shell = True)

        shutil.move(tmpfile, vepannotation)

    subprocess.check_call('export JYTHONPATH=/hgsc_software/cancer-analysis/halotron/illumina/illumina_v0.0.2/halotron/site-packages && export CLASSPATH=/hgsc_software/cancer-analysis/code/javalib/sqlite-jdbc-3.8.11.2/sqlite-jdbc-3.8.11.2.jar:/hgsc_software/cancer-analysis/code/picard-tools-1.129/picard.jar:/hgsc_software/cancer-analysis/code/picard-tools-1.129/picard-lib.jar:/hgsc_software/cancer-analysis/code/picard-tools-1.129/htsjdk-1.129.jar:/hgsc_software/cancer-analysis/code/krcgtk/krcgtk-0.01.jar:$CLASSPATH && /stornext/snfs2/can/code/jython/jython-2.7.0/bin/jython -J-Xmx15g /hgsc_software/cancer-analysis/halotron/illumina/illumina_v0.0.2/halotron/illumina/annotation/annotate_vcf_cosmic.py --reference %(reference)s %(cosmic)s %(input)s %(valstatus)s %(output)s' % {
        'reference': '/hgsc_software/cancer-analysis/resources/references/human/hg19/hg19.fa',
        'cosmic': '/hgsc_software/cancer-analysis/resources/annotation-databases/cosmic/v71/CosmicCodingMuts.cnt3.vcf',
        'input': vepannotation,
        'valstatus': '/hgsc_software/cancer-analysis/resources/dbsnp/hg19/146/dbSNP_b146_GRCh37p13.valstatus.db',
        'output': tmpfile},
        shell = True)

    logger.info("Processing %s complete", outputfpath)
    shutil.move(tmpfile, outputfpath)
    return outputfpath

def convert(opath, fpath):
    outputfpath = opath
    tmpfile = os.path.join(os.path.dirname(outputfpath), 'tmpfile.maf')
    if os.path.isfile(outputfpath):
        logger.info("Skipping conversion because %s exists", outputfpath)
        return outputfpath
    logger.info("Processing conversion of %s -> %s", fpath, outputfpath)
    nid, tid, nbar, tbar = getTNids(fpath)

    subprocess.check_call('/hgsc_software/perl/perl-5.16.2/bin/perl %(PACKAGEDIR)s/vcf2maf/vcf2maf.pl -no-annotate -input-vcf %(input)s -output-maf %(output)s -vcf-tumor-id %(tid)s -vcf-normal-id %(nid)s -tumor-id %(tbar)s -normal-id %(nbar)s -copythrough COSMIC,CENTERS,CONTEXT,DBVS' % {
        'input': fpath,
        'output': tmpfile,
        'tid': tid,
        'nid': nid,
        'tbar': tbar,
        'nbar': nbar,
        'PACKAGEDIR': PACKAGEDIR},
        shell = True)

    shutil.move(tmpfile, outputfpath)
    return outputfpath

def main(args):
    try:
        if len(args.vcfs) != len(args.callers):
            raise ValueError("vcfs and callers lengths are not the same")
        if args.tmpdir:
            if not os.path.isdir(args.tmpdir):
                os.makedirs(args.tmpdir)
        else:
            args.tmpdir = tempfile.mkdtemp()

        # generate the caller tuples
        calls = zip(args.callers, args.vcfs)
        # filter the vcf files
        filters = [(c, filter(f, c, args.tmpdir)) for c, f in calls]
        # sort the vcf files
        sorts = [(c, sort(f, args.tmpdir)) for c, f in filters]
        # v2v
        v2vs = [(c, v2v(f, args.tmpdir)) for c, f in sorts]
        # merge
        merged = merge(os.path.join(args.tmpdir, 'merged.vcf'), v2vs)
        # annotate
        annotated = annotate(merged, args.tmpdir)
        # vcf2maf
        convert(args.OUTPUTMAF, annotated)
        logger.info("Done")
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = traceback.format_exception(exc_type, exc_value, exc_traceback)
        message = 'Error in running merge.py (%s)\n\nTraceback was:\n%s\n' % (str(sys.argv), ' '.join(tb))
        notify(message)
        raise


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--vcfs', type = str, nargs = '+', help = 'vcf file path(s)')
    parser.add_argument('--callers', type = str, nargs = '+', help = 'caller keys, same len as vcfs and in same order')
    parser.add_argument('--tmpdir', type = str, help = 'location of tmp directory for processing')
    parser.add_argument('OUTPUTMAF', type = str, help = 'output file path for the merged MAF file')

    args = parser.parse_args()

    main(args)



