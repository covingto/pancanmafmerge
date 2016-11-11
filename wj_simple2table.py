

import os, os.path, sys
import json, re
import hgsc_vcf, csv

def get_allele_count(record, altnumber, samples = None):
    gt = hgsc_vcf.split_gt(record)
    if not samples:
        samples = record['SAMPLES']
    for i, g in enumerate(gt):
        if g == altnumber:
            return sum([int(s['AC'][i]) for s in samples.values()])

def new_valbatch(record):
    pos, ref, alt = extract_org_allele(record)
    valbatch = ValBatch(record['INFO']['OC'][0], pos, ref, alt)
    valbatch.add_record(record)
    return valbatch


class ValBatch(object):
    def __init__(self, oc, pos, ref, alt):
        self.oc = oc
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.records = []

    def oc_match(self, oc):
        return oc == self.oc

    def add_record(self, record):
        assert record['INFO']['OC'][0] == self.oc, "Attempt to add a record that is not part of this batch"
        self.records.append(record)

    def get_records(self):
        return self.records

    def clean(self):
        # removes records based on the clean criteria
        # if there is a single record whos sum of alts is greater than 5X the next best match
        # or whos sum allele count is greater than 2
        # then only that record is retained
        if len(self.records) < 2:
            return # we don't do anything
        varcounts = [(i, get_allele_count(r, 1)) for i, r in enumerate(self.records)]
        varcountssorted = sorted(varcounts, key=lambda x: x[1])
        if varcountssorted[-1][1] > varcountssorted[-2][1] * 5:
            max_alt_index = varcountssorted[-1][0]
            gt2_indexs = [i for i, c in varcounts if c > 2]
            self.records = [r for i, r in enumerate(self.records) if i == max_alt_index or i in gt2_indexs] # reset the list based on the max index

def get_samples(record, pat):
    return {k:v for k, v in record['SAMPLES'].items() if re.search(pat, k) is not None}

def extract_key(record):
    oc = record['INFO']['OC']
    return ';'.join([re.search(r'\[(.*?)\]', o).group(1) for o in oc])

def extract_org_allele(record):
    return parse_oc(record['INFO']['OC'][0])

def parse_oc(oc):
    g = re.search(r'\[([0-9]+)~([ACGT]+)>([ACGT]+)\]', oc)
    return g.group(1), g.group(2), g.group(3)

def main(args):
    reader = hgsc_vcf.Reader(args.INFILE)
    writer = csv.DictWriter(args.OUTFILE, delimiter = '\t', 
            fieldnames = [
                'CHROM',
                'INPOS',
                'INREF',
                'INALT',
                'VALPOS',
                'VALREF',
                'VALALT',
                'SUM_TUMOR_REF',
                'SUM_TUMOR_ALT',
                'SUM_TUMOR_DP',
                'SUM_NORMAL_REF',
                'SUM_NORMAL_ALT',
                'SUM_NORMAL_DP',
                'SUM_TUMOR_VAL_REF',
                'SUM_TUMOR_VAL_ALT',
                'SUM_TUMOR_VAL_DP',
                'SUM_NORMAL_VAL_REF',
                'SUM_NORMAL_VAL_ALT',
                'SUM_NORMAL_VAL_DP', 
                'VALKEY'])
    writer.writeheader()
    def batch(reader):
        first = reader.next()
        valbatch = new_valbatch(first)
        for record in reader:
            if not valbatch.oc_match(record['INFO']['OC'][0]):
                yield valbatch
                valbatch = new_valbatch(record)
            else:
                valbatch.add_record(record)
        yield valbatch
    for b in batch(reader):
        b.clean()
        for r in b.get_records():
            tumor_val_samples = get_samples(r, 'TUMOR_VALIDATION')
            tumor_samples = get_samples(r, r'^TUMOR$')
            normal_val_samples = get_samples(r, 'NORMAL_VALIDATION')
            normal_samples = get_samples(r, r'^NORMAL$')
            o = {
                    'CHROM': r['CHROM'],
                    'INPOS': b.pos,
                    'INREF': b.ref,
                    'INALT': b.alt,
                    'VALPOS': r['POS'],
                    'VALREF': r['REF'],
                    'VALALT': ','.join(r['ALT']),
                    'SUM_TUMOR_REF': str(get_allele_count(r, 0, tumor_samples)),
                    'SUM_TUMOR_ALT': str(get_allele_count(r, 1, tumor_samples)),
                    'SUM_TUMOR_DP': str(sum([int(s['DP'][0]) for s in tumor_samples.values()])),
                    'SUM_NORMAL_REF': str(get_allele_count(r, 0, normal_samples)),
                    'SUM_NORMAL_ALT': str(get_allele_count(r, 1, normal_samples)),
                    'SUM_NORMAL_DP': str(sum([int(s['DP'][0]) for s in normal_samples.values()])),
                    'SUM_TUMOR_VAL_REF': str(get_allele_count(r, 0, tumor_val_samples)),
                    'SUM_TUMOR_VAL_ALT': str(get_allele_count(r, 1, tumor_val_samples)),
                    'SUM_TUMOR_VAL_DP': str(sum([int(s['DP'][0]) for s in tumor_val_samples.values()])),
                    'SUM_NORMAL_VAL_REF': str(get_allele_count(r, 0, normal_val_samples)),
                    'SUM_NORMAL_VAL_ALT': str(get_allele_count(r, 1, normal_val_samples)),
                    'SUM_NORMAL_VAL_DP': str(sum([int(s['DP'][0]) for s in normal_val_samples.values()])),
                    'VALKEY': extract_key(r)
                }
            writer.writerow(o)
                
            

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('INFILE', type = argparse.FileType('r'), help = 'input vcf file')
    parser.add_argument('OUTFILE', type = argparse.FileType('w'), help = 'output vcf file')

    args = parser.parse_args()
    main(args)
