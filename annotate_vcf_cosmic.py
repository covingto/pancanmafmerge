'''
Created on Sep 9, 2015

@author: covingto
'''

import os
import os.path
import sys
#TODO: refactor to use hgsc_vcf library
import hgsc_vcf
import re
import logging
import htsjdk.samtools.reference.IndexedFastaSequenceFile as IndexedFastaSequenceFile
import org.bcm.hgsc.utils.SynchronousIndexedFastaReader as SFastaReader
import java.io.File as File
import java.lang.String as String
from java.lang import Class
from java.sql  import DriverManager, SQLException

JDBC_DRIVER = "org.sqlite.JDBC"

def getConnection(jdbc_url, driverName):
    """
        Given the name of a JDBC driver class and the url to be used 
        to connect to a database, attempt to obtain a connection to 
        the database.
    """
    try:
        Class.forName(driverName).newInstance()
    except Exception, msg:
        print msg
        sys.exit(-1)

    try:
        dbConn = DriverManager.getConnection(jdbc_url)
    except SQLException, msg:
        print msg
        sys.exit(-1)

    return dbConn

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

re_csqformat = re.compile(r'.*?Format: (.*)')
re_capstring = re.compile(r'[A-Z]+')
##
# extract csq format information from a vcf info description


def get_csq_format(csq_descriptor):
    _re_csq = re.search(re_csqformat, csq_descriptor)
    if _re_csq is None:
        raise ValueError("Format search failed")
    return _re_csq.group(1).split('|')

##
# codon types
SIMPLE = 0B0001
CODESHIFTING = 0B0010
INS = 0B0110
DEL = 0B1010

##
# returns the offset for the codon range (if one exists) for this csq codon string
#
# given a string, this locates the offset of the capital letter and the next offset
# for the end of the string, thereby defining the codon
# The string gTt/gCt would return -1 and +1 indicating that the codon begins one
# base before the allele and ends one base after it
# The following offsets may be observed for coding snv's
# Aaa/Taa: 0, +2
# aAa/aTa: -1, +1
# aaA/aaT: -2, 0
#
# The following is used for deletions
# ggAGGg/ggg: -2, +3
#
# The following is used for insertions
# a|cg/aCGCcg: 0, +2
# in this case the bar represents the location of the insertion (between the a and c)
# the coordinates in this case match the VCF format with the insertion involving the first base of the insertion
# note that the insertion will only ever involve a range of 3 bases (inclusively)
#
# complex substitutions (sometimes annotated as del/ins) are handled be returning the
# offsets relative to the start of the allele as in the following:
# cAGGTGCTg/cGGGTTCTg: -1, +7
# 101234567          : key


def get_cqs_codon_region(codon_s):
    _c1, _c2 = codon_s.split('/')
    # codon deletion and insertion events
    if _c1 == '-':
        return -2, 4, 0, SIMPLE + INS
    elif _c2 == '-':
        return 0, 0, len(_c1), SIMPLE + DEL
    _c1_splitter = re.findall(re_capstring, _c1)
    # first test _c1 since that is where most of the action is
    if len(_c1_splitter) > 1:
        raise ValueError(
            "Interspersed capitalization for codon1 in %s" %
            codon_s)
    elif len(_c1_splitter) < 1:  # move to _c2, implies INS
        _c2_splitter = re.findall(re_capstring, _c2)
        if len(_c2_splitter) != 1:
            raise ValueError(
                "Codon1 and Codon2 have no capital letters for splitting or interspersed capitalization in %s" %
                codon_s)
        _left, _right = _c2.split(_c2_splitter[0])
        return abs(len(_left) - 1) * -1, len(_right) + 1, 0, INS
    else:  # implies either DEL or SIMPLE
        _left, _right = _c1.split(_c1_splitter[0])
        if len(_c1) == len(_c2):
            _type = SIMPLE
        else:
            _type = DEL
        return len(_left) * -1, len(_right), len(_c1_splitter[0]), _type

##
# returns the genomic range of a codon
#
# given an allele start position and the codon alterations for that allele, this
# function returns the genomic range of bases that comprise the codon(s) altered
# by the allele
#
# note that the positions listed below are in VCF coordinates
# while the codon resolver is in 1 base (the location of the alternation)
# that means that deleted alleles start at the base that was removed
# insertions happen at the preceeding base
# 100    A    T       , Aaa/Taa:      100, 102
# 101    A    T       , aAa/aTa:      100, 102
# 102    A    T       , aaA/aaT:      100, 102
# 99    gGTA  g       , GTA/-:        100, 102
# 99    g     gGTA    , -/GTA:         97, 102  # codon boundary case; return codons before and after
# 100   g     gAAA    , gta/gAAAta:   100, 102
# 100   gTAG  g       , gTAGgg/ggg:   100, 105


def get_codon_genomic_range(position_i, codon_s, pos_strand_b=True):
    _left, _right, _bases_altered, _codon_type_i = get_cqs_codon_region(
        codon_s)

    if not pos_strand_b:
        if not (_codon_type_i & CODESHIFTING == CODESHIFTING) or not (
                _codon_type_i & SIMPLE == SIMPLE):
            _left, _right = _right * -1, _left * -1
            if (_codon_type_i & INS == INS):
                # compensate for inclusivity and zero basing
                _left += 2
                _right += 2
    _right += _bases_altered - 1
    return position_i + _left, position_i + _right



##
# Conatiner for VCF files
#
# allows rapid searching and manipulation of VCF files in memory
class VCFContainer:

    def __init__(self, reader, buffer):
        self.reader = reader
        self.buffer = buffer
        self.records = {}
        self.__load()
        self.lastindex = 0
        self.lastpos = 0

    def __load(self):
        logger.info("Reading VCF file into memory")
        for record in self.reader:
            if record['CHROM'] not in self.records:
                self.records[record['CHROM']] = []
            self.records[record['CHROM']].append(record)

    ##
    # return a list of records within the indicated buffer limit
    def intersect(self, record):
        if not record['CHROM'] in self.records:
            return []
        _record_set = self.records[record['CHROM']]
        _pos = record['POS'] - self.buffer
        _end = record['POS'] + len(record['REF']) + self.buffer

        if _pos > self.lastpos:
            self.lastindex = 0
        self.lastpos = _pos

        _result = []
        for _i in range(self.lastindex, len(_record_set)):
            _test_record = _record_set[_i]
            if ((_test_record['POS'] + len(_test_record['REF']))
                    >= _pos) and (_test_record['POS'] <= _end):
                _result.append(_test_record)
            elif _test_record['POS'] > _end:
                break  # exceeded the end, stop
            self.lastpos = _test_record['POS']
            self.lastindex = _i

        return _result

MATCH_PRIORITY = {'SITE': 0, 'CODON': 1, 'BUFFER': 2, 'NONE':3}
##
# match up cosmic records and format for printing
#
# matches are records that are within the buffer of the indicated site for the primary record
# csq are the csq info fields (may be multiple) corresponding to this record and it's alts
# the goal is to convert each of the matches to an appropriate string.  The function will return
# the most severt csq match for each of the matches in the order SITE,
# CODON, BUFFER
def generate_cosmic_info(matches, csq_l, record):
    if len(matches) == 0:
        return ['NONE']
    _result = []
    for _m in sorted(matches, key = lambda _m: MATCH_PRIORITY[get_match_type(_m, csq_l, record)]):
        _result.append('|'.join([get_match_type(_m,
                                                csq_l,
                                                record),
                                 _m['INFO']['AA'][0],
                                 _m['INFO']['CDS'][0],
                                 str(_m['INFO']['CNT'][0])]))
    return _result

##
# return the match type for this match and csq list
#
# match is a single record form a cosmic vcf and csq is a list of csq fields from vep
# the default of no other match is found is to report
def get_match_type(match, csq_l, record):
    _sites = [get_csq_cdna_site(_csq) for _csq in csq_l]
    if match['INFO']['CDS'][0] in _sites:
        return 'SITE'
    _codon_ranges = [
        get_codon_genomic_range(
            record['POS'],
            _csq['Codons'],
            _csq['STRAND'] == 1) for _csq in csq_l if _csq['Codons'] != '']
    for _start, _end in _codon_ranges:
        if match['POS'] < _end and (match['POS'] + len(match['REF']) > _start):
            return 'CODON'
    return 'BUFFER'


def get_csq_cdna_site(csq):
    try:
        return csq['HGVSc'].split(':')[1]
    except:
        return ''

def add_info_to_reader(reader, header):
    reader.header.add_header(header)

def add_command_to_reader(reader, command):
    reader.header.add_header(command)

def get_val_status(rsid, con):
    resultSet = con.executeQuery("select rsid, valstatus from dbsnpvalstat where rsid = \"%s\"" % rsid)
    while resultSet.next():
        yield resultSet.getString("valstatus")

def generate_valstatus_info(_existing_ids, _valdata):
    dbCon = getConnection("jdbc:sqlite:%s"  % _valdata, JDBC_DRIVER)
    con = dbCon.createStatement()
    try:
        _vals = [v for _id in _existing_ids for v in get_val_status(_id, con)]
        _vals_set = set(_vals)
        if len(_vals_set) > 0:
            return "|".join(_vals_set)
        else:
            return '.'
    finally:
        con.close()
        dbCon.close()

def main(args):
    vcf_reader = hgsc_vcf.Reader(args.INPUTVCF)
    vcf_container_cosmic = VCFContainer(
        hgsc_vcf.Reader(
            args.COSMICVCF),
        args.buffer)
    # read in the dbsnp data
    # connect to the reference file
    ifasta = IndexedFastaSequenceFile(File(args.reference))
    add_command_to_reader(vcf_reader, '##COMMAND=<ID=annotate_vcf_cosmic.py,Params="%s">' % " ".join(sys.argv))
    # add the COSMIC header info
    add_info_to_reader(vcf_reader, '##INFO=<ID=COSMIC,Number=.,Type=String,Description="' + 
        'COSMIC info, can be one of NONE, BUFFER, CODON, SITE.  ' +
        'All but NONE are accompanied by AA|CDS|CNT BUFFER indicates the COSMIC site is within %(buffer)sbp of the position.  example: ' +
        'SITE|p.P228fs*227|c.682_683insT|3 or NONE.  VCF file used was %(cosmicvcf)s.">\n' % {'buffer': str(args.buffer), 'cosmicvcf': args.COSMICVCF})

    # add the context
    add_info_to_reader(vcf_reader, "##INFO=<ID=CONTEXT,Number=1,Type=String,Description=\"Base context around variant. [POS - 5, POS + len(REF) + 4]\">\n")
    # add the validatio status info
    add_info_to_reader(vcf_reader, "##INFO=<ID=DBVS,Number=1,Type=String,Description=\"dbSNP validation status, | separated\">\n")
    # get the format for the vep annotations
    _vep_format = get_csq_format([h for h in vcf_reader.header.get_headers('INFO', 'CSQ')][0].fields['Description'])

    vcf_writer = hgsc_vcf.Writer(args.OUTPUTVCF, vcf_reader.header)
    vcf_writer.write_header()
    for record in vcf_reader:
        try:
            ## check that the position is annotated with CSQ, if not then this is a write through
            if 'CSQ' in record['INFO']:
                # matches are intersecting hits in the VCF
                
                _matches = vcf_container_cosmic.intersect(record)
                _csq_l = [dict(zip(_vep_format, _csq.split('|')))
                          for _csq in record['INFO'].get('CSQ')]
                _info = generate_cosmic_info(_matches, _csq_l, record)
                record['INFO']['COSMIC'] = _info
                # extract the dbsnp validation rsids
                _existing_ids = [_id for _csq in _csq_l for _id in _csq['Existing_variation'].split('&')]
                record['INFO']['DBVS'] = [generate_valstatus_info(_existing_ids, args.DBSNPVAL)]
            record['INFO']['CONTEXT'] = [str(
                    String(
                        ifasta.getSubsequenceAt(record['CHROM'], record['POS'] - 5, record['POS'] + len(record['REF']) + 4).getBases()
                        ))]
            
        except:
            logger.exception("Error in record modification")
        vcf_writer.write_record(record)


##
# returns the amino acid residue number from ...


def get_csq_hgvsp_aa_number(csq_hgvsp_s):
    pass

if __name__ == '__main__':
    # set up the logger
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    logger.addHandler(ch)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--reference',
            type = str,
            help = 'Reference file for context annotation'
            )
    parser.add_argument(
        '--buffer',
        default=10,
        type=int,
        help='Buffer size for matching COSMIC records, defaults to 10bp')
    parser.add_argument(
        'COSMICVCF',
        type=argparse.FileType('r'),
        help='COSMIC vcf for annotation')
    parser.add_argument(
        'INPUTVCF',
        type=argparse.FileType('r'),
        help='Input vcf to annotate')
    parser.add_argument('DBSNPVAL',
                        type = str,
                        help = "dbSNP validation data in the form rsid [tab] valstatus")
    parser.add_argument(
        'OUTPUTVCF',
        type=argparse.FileType('w'),
        nargs='?',
        default=sys.stdout,
        help='Output file, defaults to stdout')

    args = parser.parse_args()

    main(args)
