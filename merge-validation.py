

import os, os.path, sys
import collections
import gzip
import csv, json

def update_validation(maf_record, val_dict):
    # maf_record is coming from the input maf file
    # if we find a match to something in the val_dict, we pull that info and also mark the val_dict
    # record as being hit, this lets us check how many validation inputs were missed
    # we will write the val_dict back out as json in case we need to check it later

    # we will set this up so that the val_dict is first indexed by subject and then by valkey which 
    # then points to the validation info (Validation_Status, Mutation_Status, Validation_Method

    # we modify the maf_record in place, so no return here

    if maf_record.get('Tumor_Sample_Barcode')[:12] not in val_dict:
        return
    valkey = make_val_key(maf_record)
    valinfo = val_dict[maf_record.get('Tumor_Sample_Barcode')[:12]].get(valkey, None)
    if valinfo is None:
        return
    else:
        maf_record.update(valinfo)
        if valinfo['orgbarcode'] != maf_record.get('Tumor_Sample_Barcode'):
            maf_record['Validation_Status'] = maf_record['Validation_Status'] + '*'
        valinfo['hit'] = True
        return

def make_val_key(maf_record):
    valkey = '|'.join([maf_record['Chromosome'], maf_record['Start_Position'], maf_record['End_Position'], maf_record['Reference_Allele'], maf_record['Tumor_Seq_Allele2']])
    return valkey

def build_val_dict(fpath):
    # the input is a validation maf file
    result = collections.defaultdict(dict)
    with open(fpath, 'r') as fi:
        reader = csv.DictReader(fi, delimiter = '\t')
        for r in reader:
            result[r['Tumor_Sample_Barcode'][:12]][make_val_key(r)] = {k:v for k, v in r.items() if k in ('Validation_Status', 'Mutation_Status', 'Validation_Method')}
            result[r['Tumor_Sample_Barcode'][:12]][make_val_key(r)]['orgbarcode'] = r['Tumor_Sample_Barcode']
    return result

def validate_maf(maf_path, val_path, out_path, val_used_path):
    print "building validation dict"
    val_dict = build_val_dict(val_path)
    print "processing maf file"
    with gzip.open(maf_path, 'rb') as fi, gzip.open(out_path, 'wb') as fo:
        reader = csv.DictReader(fi, delimiter = '\t')
        writer = csv.DictWriter(fo, delimiter = '\t', fieldnames = reader.fieldnames, extrasaction = 'ignore')
        writer.writeheader()
        for i, r in enumerate(reader):
            update_validation(r, val_dict)
            writer.writerow(r)
            if i % 10000 == 0:
                print "processed %s records" % str(i)
    with open(val_used_path, 'w') as fo:
        json.dump(val_dict, fo)
    used = total = 0
    for v in val_dict.values(): # these are the dicts for the subjects
        for vv in v.values(): # these are the dicts for the validations
            total += 1
            if vv.get('hit', False):
                used += 1
    print "Used %s of %s total validations" % (str(used), str(total))

def test():
    vkeytest = make_val_key({'Chromosome': '1', 'Start_Position': '100', 'End_Position': '101', 'Reference_Allele': 'A', 'Tumor_Seq_Allele2': 'C'})
    assert vkeytest == '1|100|101|A|C', "didn't make correct val_key"
    val_dict_test = {'TCGA-AA-1234': {vkeytest: {'orgbarcode': 'TCGA-AA-1234-05', 'Validation_Status': 'Valid', 'Mutation_Status': 'Somatic', 'Validation_Method': 'Sanger'}}}
    record1test = {'Tumor_Sample_Barcode': 'TCGA-AA-1234-05', 'Chromosome': '1', 'Start_Position': '100', 'End_Position': '101', 'Reference_Allele': 'A', 'Tumor_Seq_Allele2': 'C',
            'Validation_Status': 'unknown', 'Mutation_Status': 'Somatic', 'Validation_Method': '.'}
    record2test = {'Tumor_Sample_Barcode': 'TCGA-AA-1234-05', 'Chromosome': '1', 'Start_Position': '110', 'End_Position': '111', 'Reference_Allele': 'A', 'Tumor_Seq_Allele2': 'C',
            'Validation_Status': 'unknown', 'Mutation_Status': 'Somatic', 'Validation_Method': '.'}
    record3test = {'Tumor_Sample_Barcode': 'TCGA-AA-1234-06', 'Chromosome': '1', 'Start_Position': '100', 'End_Position': '101', 'Reference_Allele': 'A', 'Tumor_Seq_Allele2': 'C',
            'Validation_Status': 'unknown', 'Mutation_Status': 'Somatic', 'Validation_Method': '.'}
    update_validation(record1test, val_dict_test)
    update_validation(record2test, val_dict_test)
    update_validation(record3test, val_dict_test)
    assert record1test['Validation_Status'] == 'Valid', "didn't update record1 to be Valid"
    assert record2test['Validation_Status'] == 'unknown', "changed the value of an unmatched record"
    assert val_dict_test['TCGA-AA-1234'][vkeytest]['hit'] == True, "did not record the hit"
    assert record3test['Validation_Status'] == 'Valid*', "didn't get that the validation status was a partial hit"
    print "passed all tests, hope this works :)"

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--test', action = 'store_true', help = 'run a quick test of the functions')
    parser.add_argument('--mafpath', type = str, help = 'path to the gzip maf file')
    parser.add_argument('--valpath', type = str, help = 'path to the validation maf file')
    parser.add_argument('--outpath', type = str, help = 'path to the output validated gzip maf file')
    parser.add_argument('--valusedpath', type = str, help = 'path to the validation dict result file')

    args = parser.parse_args()
    if args.test:
        test()
    else:
        validate_maf(args.mafpath, args.valpath, args.outpath, args.valusedpath)
