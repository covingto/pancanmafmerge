##
# @author Kyle R. Covington; kylecovington1@gmail.com
# 
# @purpose To dispatch merge.py for all samples in a study


import os, os.path, sys
import glob
import subprocess
import logging

logger = logging.getLogger('run-batch')
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler(stream=sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

logger.addHandler(ch)

PACKAGE = os.path.dirname(os.path.abspath(__file__))

##
# dispatch a job to the infrastructure
#
# @param dpath - the path to a vcf directory.  The directory must contain files that glob match the 7 globs for pancan vcf files.
# @returns None
# @effects: if the output file path (resultpath/`basename dpath`/merged.maf) exists then nothing happens, otherwise required dirs are created and merge.py is executed
def dispatch(dpath, resultpath):
    logger.info("Processing %s", dpath)
    outdir = os.path.join(resultpath, os.path.basename(dpath))
    output = os.path.join(outdir, 'merged.maf')
    tmpdir = os.path.join(outdir, 'tmp')
    if os.path.isfile(output):
        logger.info("%s has already been generated, remove if you want to run again", output)
        return
    else:
        for d in (outdir, tmpdir):
            if not os.path.isdir(d):
                os.makedirs(d)
    # get the callers
    vcfs = []
    callers = []
    for gpat in ('radia', 'muse', 'mutect', 'pindel', 'SomaticSniper', 'varscan.indel', 'varscan.snp'):
        gfiles = glob.glob(os.path.join(dpath, '*%s*' % gpat))
        if len(gfiles) != 1:
            logger.error("did not find 1 file for %s: %s", gpat, gfiles)
            return
        vcfs.append(gfiles[0])
        if gpat == 'varscan.indel':
            callers.append('VARSCANI')
        elif gpat == 'varscan.snp':
            callers.append('VARSCANS')
        else:
            callers.append(gpat.upper())
    
    CMD = 'echo "python %(PACKAGE)s/merge.py --vcfs %(vcfs)s --callers %(callers)s --tmpdir %(tmpdir)s %(output)s" | msub -V -A proj-dwpancan -q analysis -N %(sample)s -d $PWD -j oe -o %(output)s.log -l nodes=1:ppn=4,mem=32g' % {
        'output': output,
        'sample': os.path.basename(dpath),
        'tmpdir': tmpdir,
        'callers': ' '.join(callers),
        'vcfs': ' '.join(vcfs),
        'PACKAGE': PACKAGE}
    logger.info(CMD)
    subprocess.check_call(CMD, shell = True)
    # logger.info("Dispatched job for %s", dpath)

def main(args):
    for dpath in args.indir:
        if os.path.isdir(dpath):
            dispatch(dpath, args.resultdir)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--resultdir', type = str, help = 'directory where results will be stored')
    parser.add_argument('indir', nargs = '+', type = str, help = 'pair vcf directories')

    args = parser.parse_args()

    main(args)
