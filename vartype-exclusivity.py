'''
Created on Feb 4, 2016

@author: covingto
'''


import hgsc_vcf
import os
import sys
import re, json
import base64, hashlib
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.DEBUG)

def safe_div(n, d):
    try:
        return float(n) / float(d)
    except ZeroDivisionError:
        return 0

def merge_map(x, y):
    newmap = {}
    newmap.update(x)
    for k, v in y.items():
        if k not in newmap:
            newmap[k] = v
        elif isinstance(v, list) and isinstance(newmap[k], list):
            newmap[k] += v
        elif isinstance(v, dict) and isinstance(newmap[k], dict):
            newmap[k] = merge_map(newmap[k], v)
        else:
            raise Exception("maps contain non-conformant: %s and %s", x, y)
    return newmap


##
# Used for making the final sample model map.
#
# Returns a dict of nearly the same structure as the @param _model
# but with regex terms replaced with matches to sample names
# in the @param samples list.  The regex terms are matched against
# roles in the @param _map.
def make_sample_direct_map(_map, _model):
    # _model is a dict, if its value is a dict, recurse
    if isinstance(_model, dict):
        tmp = {k:make_sample_direct_map(_map, v) for k, v in _model.items()}
        return {k:v for k, v in tmp.items() if v}
    elif isinstance(_model, list):
        matches = []
        for regex in _model:
            matches += [k for k, v in _map.items() if regex == k]
        return matches

##
# Used for making the final sample model map.
#
# Returns a dict of nearly the same structure as the @param _model
# but with regex terms replaced with matches to sample names
# in the @param samples list.  The regex terms are matched against
# roles in the @param _map.
def make_sample_map(_map, _model):
    # _model is a dict, if its value is a dict, recurse
    if isinstance(_model, dict):
        tmp = {k:make_sample_map(_map, v) for k, v in _model.items()}
        return {k:v for k, v in tmp.items() if v}
    elif isinstance(_model, list):
        matches = []
        for regex in _model:
            newmatches = []
            for k, v in _map.items():
                try:
                    if re.match(regex + "$", v):
                        newmatches.append(k)
                except:
                    logger.warning("regex problem with %s and %s", regex, v)
            matches += newmatches
        return matches


## 
# generate sample pairs from a model_map
#
# effectively recurses and selects the most basic elements from the list
def get_samples_from_model_map(model_map):
    if isinstance(model_map, dict):
        for k, v in model_map.items():
            for s in get_samples_from_model_map(v):
                if isinstance(s, basestring):
                    yield s, k 
                else:
                    yield s
    elif isinstance(model_map, list):
        for s in model_map:
            yield s 

##
# build the affinity function
#
# the config looks like:
#
# {'sample': {
#    'map': {
    #    2: "r['MRQ'][altindex] > 1 and r['AQC'][altindex] > 2 and float(r['AC'][altindex])/float(r['DP']) > 0.04 and r['FC'][altindex] > 1 and r['RC'][altindex] > 1",
    #    1: "float(r['AQC'][altindex])/float(r['DP']) > 0.01",
    #    0: "r['AC'][altindex] < 1 and r['DP'] > 10",
    #    -1: "True"
#    }, 
#    'order': [2,1,0,-1]
#    }
# }
def build_sample_affinity_function(_config):
    _sample_config = _config.get('sample')
    def _sample_affinity(r, refindex, altindex):
        for k in _sample_config['order']:
            try:
                if eval(_sample_config['map'][k]):
                    return int(k) # must convert back to an int since ints are not valid keys in json
            except Exception as inst:
                print r
                raise ValueError("Error in evaluation of %s: %s", _sample_config['map'][k], str(inst))
    return _sample_affinity

def build_record_mod_function(_config):
    def _record_mod(r):
        mods = []
        for k, v in _config['record'].items():
            try:
                if eval(v):
                    mods.append(k)
            except Exception as inst:
                print r
                raise ValueError("Error in evaluating %s: %s", v, str(inst))
        return ''.join(mods)
    return _record_mod

def get_terminal_samples(d):
    for k, v in d.items():
        if isinstance(v, dict):
            for t in get_terminal_samples(v):
                yield t
        else:
            for t in v:
                yield t

def main(args):
    if os.path.isfile(args.MODEL):
        with open(args.MODEL, 'r') as fi:
            _model = json.load(fi)
    else:
        _model = json.loads(args.MODEL)

    if os.path.isfile(args.TYPEMAP):
        with open(args.TYPEMAP, 'r') as fi:
            _map = json.load(fi)
    else:
        _map = json.loads(args.TYPEMAP)
    logger.info("Input map: %s", _map)
    _config = json.load(args.CONFIG)

    vcf = hgsc_vcf.Reader(args.INPUT)
    
    # reduce the map to just the samples used in this study
    _map = {k:v for k, v in _map.items() if k in vcf.header.samples}
    logger.info("Revised map: %s", _map)
    logger.info("Direct Mapping")
    sample_model_direct_map = make_sample_direct_map(_map, _model)
    logger.info(sample_model_direct_map)
    sample_type_direct_mapping = set([v for v in get_samples_from_model_map(sample_model_direct_map)])
    logger.info(sample_type_direct_mapping)
    sample_model_regex_map = make_sample_map({k:v for k, v in _map.items() if k not in [s for s, m in sample_type_direct_mapping]}, _model)
    logger.info(sample_model_regex_map)
    # check that there are no duplicates
    sample_type_regex_mapping = set([v for v in get_samples_from_model_map(sample_model_regex_map)])
    logger.info("Regex Mapping")
    logger.info(sample_type_regex_mapping)
    # take the mappings that are direct first, then the ones that are satisfied via regex
    direct_map_samples = [s for s, k in sample_type_direct_mapping]
    regex_only_mappings = [(s, k) for s, k in sample_type_regex_mapping if s not in direct_map_samples]
    sample_type_mapping = list(sample_type_direct_mapping) + list(set(regex_only_mappings))
    sample_model_map = merge_map(sample_model_direct_map, sample_model_regex_map)
    logger.info("Final Mapping")
    logger.info(sample_model_map)
    logger.info(sample_type_mapping)
    _used_samples = [s for s in sample_type_mapping]
    if len(_used_samples) != len(set(_used_samples)):
        raise ValueError("You have samples that are listed in two branches, indicates a role collision: %s", _used_samples)
    if len(_used_samples) != len(vcf.header.samples):
        logger.error("Used samples: %s", _used_samples)
        logger.error("Header samples: %s", vcf.header.samples)
        _us = [v[0] for v in _used_samples]
        for s in vcf.header.samples:
            if s not in _us:
                logger.error("No match for %s",s)
        raise ValueError("You have not mapped all samples in the vcf, check that the roles for each sample in the TYPEMAP are used in the MODELMAP")
        
    
    # build the _sample_affinity_f
    _sample_affinity_f = build_sample_affinity_function(_config)
    _record_f = build_record_mod_function(_config)
    
    writer = hgsc_vcf.Writer(args.OUTPUT, vcf.header)
    writer.header.add_header('##INFO=<ID=VTE,Number=1,Type=String,Description="Variant type exclusivity based on the input models and mappings.  If multiple types are detected they are separated by |.">')
    writer.header.add_header('##FORMAT=<ID=VTES,Number=1,Type=Integre,Description="Variant type exclusivity score">')
    writer.header.add_header('##COMMAND=<ID=varytpe-exclusivity,ARGS="%s">' % re.escape(' '.join(sys.argv)))
    header_samples = {s['ID']: s for s in writer.header.get_headers('SAMPLE')}
    for s, t in sample_type_mapping:
        if s in header_samples:
            header_samples[s]['TYPE'] = t
            header_samples[s]['ROLE'] = _map[s]
        else:
            writer.header.add_header('##SAMPLE=<ID=%s,TYPE=%s,ROLE=%s>' % (s, t, _map[s]))
            header_samples = {s['ID']: s for s in writer.header.get_headers('SAMPLE')}
    if args.subject:
        subject = args.subject
    else:
        subject = base64.urlsafe_b64encode(hashlib.md5("".join(sorted(header_samples.keys()))).digest())
    if len([h for h in writer.header.get_headers('SUBJECT')]) < 1:
        writer.header.add_header('##SUBJECT=<ID="%s">' % subject.replace('=', ''))
    writer.write_header()
    for record in vcf:
        for b_record in hgsc_vcf.select_allele(record, lambda x: [hgsc_vcf.best_alt_index(x)]):
            b_record['INFO']['VTE'] = [vartype_exclusivity(sample_model_map, b_record, _sample_affinity_f) + _record_f(b_record)]
            writer.write_record(b_record)


##
# Calculate vartype exclusivity.
#
# The goal is to identify the point highest up the tree to which the
# variant has affinity.  Then the variants are "exclusive" to that
# type (or it's children).
#
# A model can look like this:
# {"_S_": {
# "_O_": { # this is the organ environment from which the test samples come
# "_TAN_": ["tan.*", "TAN.*"], # maps to all samples with a role of tan or TAN
# "_T_": ["tumor.*", "Tumor.*", "T", "M", "metastasis.*"], # maps to all samples with a role that decents from a tumor
#               },
# "_C_": { # this matches all non-affected orgain system samples
#                "_B_": ["B", "Blood.*", "blood.*"]
#               },
#        }
# }
#
# A map can look like this:
# {"S1": "TAN1",
#  "S1-R": "TAN-rna",
#  "S1-B": "blood",
#  "Tumor": "tumor",
#  "tumor-hd": "tumor-high-depth"
# }
#
# Note that all of the available mappings in the model are not used, you may make
# the model as complicated or as simple as you like, but all roles indicated in the
# map must match a role in the model.
#
# The affinity function is generated from the config and
# will return an int of -1, 0, 1, 2.
# <ul>
#    <li>-1: There is not sufficient data to support or refute affinity</li>
#    <li>0: There is sufficient data to establish affinity and this sample does not
#            have affinity for this allele</li>
#    <li>1: There is some indication that this sample has affinity for this
#            allele but it is weak and depends on a related sample also showing
#            affinity</li>
#    <li>2: There is evidence for the allele in this sample</li>
# </ul>
#
# Example:
#
#    COSMIC    S1    S1-R    S1-B    Tumor    tumor-hd    |    VTE    Explaination
#    NONE      0     0       0       2        2           |    _T_    Easy case, we only see the variant data in the _T_ derived samples, no modifiers
#    KRASG12D  0     0       0       2        2           |    _T_!   Easy case (above) ! modifier because of cosmic
#    NONE      2     2       2       2        2           |    _S_    Ease case, all samples have evidence so most ancestral is the _S_
#    NONE      -1    -1      -1      2        2           |    _T_*   Harder case, not enough evidence to make a call in more ancestral samples, _T_ anchors but modifier * indicates that we can't exclude other genetics
#    NONE      0     0       0       1        1           |    _T_-   Variation only supported by uncertain data modifier (-)
#    NONE      0     0       0       2        0           |    _T_^   Variation exclusive to a type but with conflicting modifier (^)
#
#
# @param affinity - a function that indicates if this variant has an
#                    affinity for a sample, returns -1, 0, 1, 2
def vartype_exclusivity(_map, record, _sample_affinity_f):
    samples = record['SAMPLES']
    if '.' in record['ALT']:
        return "NONE"
    
    record['FORMAT'].append('VTES')
    gt = hgsc_vcf.split_gt(record)
    refindex = hgsc_vcf.ref_index(record, gt)
    altindex = hgsc_vcf.best_alt_index(record, gt)
    # generate an affinity map for each sample
    sample_affinity = {k: _sample_affinity_f(v, refindex, altindex) for k, v in samples.items()}
    for k, v in samples.items():
        v['VTES'] = [str(sample_affinity.get(k))]
        
    # now reorganize the data to match a sample set
    # e.g. _P_:{_O_:{_TAN_:[0], _T_:[1,2]}, _C_:{_B_:[0]}} 
    # the _map will look like 
    #     _P_:{_O_:{_TAN_:[S1, S1-R], _T_:[Tumor, tumor-hd]}, _C_:{_B_:[S1-B]}}
    p = build_set(sample_affinity, _map)
    a, mod = [affinity(k, v) for k, v in p.items()][0]
    if a == [] or a is None:
        a = ['NONE']
    affin = '|'.join(a) + mod_symbol(mod)

    return affin

def build_set(sample_affinity, _map):
    if isinstance(_map, dict):
        return {k: build_set(sample_affinity, v) for k, v in _map.items()}
    elif isinstance(_map, list):
        return [sample_affinity.get(k, -1) for k in _map]

HQ = 0B0001 # the set contains at least one high quality ALT sample
RR = 0B0010 # the set contains at least one high quality RR sample
LQ = 0B0100 # the set contains at least one low quality ALT sample
AMB = 0B1000 # the existance of the variant can not be entirely excluded because there is not sufficient data at this loci


##
# returns an affinity and modifier
#
# the affinity returned is a list, which can be concatenated 
# this list component comes into play when there is a trifurcating
# branch, should a node have more than two branches and there is no
# evidence in all branches then a list with more than one entry is returned
# this may later be collapsed should higher level nodes contain evidence 
# across their branches
def affinity(k, v):
    if isinstance(v, dict):
        matches = [] # the matching branches
        m_mods = 0B0000 # mods along the matching branch
        nm_mods = 0B0000 # mods along any unmatching branch
        miss_match = False
        for _k, _v in v.items():
            affin, mod = affinity(_k, _v)
            if affin is not None:
                matches += affin
                m_mods = m_mods | mod
            else:
                miss_match = True
                nm_mods = nm_mods | mod
        # allmatch = len([m for m in v.keys() if m in matches]) == len(v.keys())
        if not miss_match:
            # remove the AMB mod from m_mods, this is beacuse we match all branches at least once
            m_mods &= ~AMB
            return [k] if len(v.items()) > 1 else matches, m_mods # | (nm_mods & AMB) # mark the set as ambiguous
        else:
            return matches if matches else None, m_mods | (nm_mods & AMB) # mark the set as ambiguous
    elif isinstance(v, list):
        mod = 0B0000
        if 2 in v:
            mod = mod | HQ
        if 1 in v:
            mod = mod | LQ
        if 0 in v:
            mod = mod | RR
        max_v = max(v)
        if max_v < 0:
            mod = mod | AMB
        if (mod & HQ) or (mod & LQ):
            return [k], mod
        else:
            return None, mod
    else:
        raise ValueError("Value v [%s] is not a list or a dict, its a %s", str(v), type(v))

def mod_symbol(mod):
    mods = []
    if not mod & HQ:
        mods.append('-')
    if mod & AMB:
        mods.append('*')
    if (mod & HQ or mod & LQ) and (mod & RR):
        mods.append('^')
    return ''.join(mods)

#[''.join(a)+mod_symbol(m) for a, m in [affinity(k, v) for k, v in s.items()]]

if __name__ == '__main__':
    mainlogger = logging.getLogger()
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    mainlogger.addHandler(ch)

    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--subject', type = str, help = 'subject name, will default to the base64, url-safe, md5 sum of the concatenation of all samples (after sorting) and is therefore deterministic if not very pretty')
    parser.add_argument(
        'MODEL', type=str, help="JSON or a .json file describing the genetic model.  The model should be described as a set of nested json objects (dict) with each key representing either a role. e.g. {\"P\":{\"B\":{}, \"T\":{}}}")
    parser.add_argument(
        'TYPEMAP', type=str, help="JSON or a .json file describing the mapping between samples and roles used in the model. e.g. {\"Sample1\": \"T\", \"Sample2\":\"B\"}")
    parser.add_argument(
        'CONFIG', type=argparse.FileType('r'), help='config file for settings')
    parser.add_argument(
        'INPUT', type=argparse.FileType('r'), help="input VCF file")
    parser.add_argument(
        'OUTPUT', type=argparse.FileType('w'), help="output VCF file")

    args = parser.parse_args()

    main(args)
