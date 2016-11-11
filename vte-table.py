

import os, os.path, sys
import logging
import hgsc_vcf
import re, csv
import bz2

logger = logging.getLogger('vte-table')
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.WARN)

re_csqformat = re.compile(r'.*?Format: (.*)')
re_capstring = re.compile(r'[A-Z]+')

STANDARD_FIELDS = ['CHROM', 'POS', 'REF', 'ALT', 'SUBJECT', 'SAMPLE', 'ROLE', 'TYPE', 'SCORE', 'VTE']

##
# pulls data from the SAMPLE info tags "##SAMPLE=<ID...>" and SUBJECT tag "##SUBJECT=<ID=...>
#
# these fields are processed and added to the output, subject is added as a primary key
# while the SAMPLE information is recoded as a map and stored for access by the 
# @ref process_vte_sample function.
# @param header - a @ref VCFHeader object
# @param output - a dict to be modified in place
def process_vte_header(header, output):
    sample_headers = [h for h in header.get_headers('SAMPLE')]
    subject_headers = [h for h in header.get_headers('SUBJECT')]
    if len(subject_headers) != 1:
        raise ValueError("There were %s SUBJECT tags, exactly 1 is required", str(len(subject_headers)))
    if len(sample_headers) < 1:
        raise ValueError("There were no SAMPLE tags")
    sample_map = {h['ID']:(h['ROLE'].strip('"'), h['TYPE'].strip('"')) for h in sample_headers}
    output['SUBJECT'] = subject_headers[0]['ID'].strip('"')
    output['_sample_map'] = sample_map

def process_vte_record(record, output):
    output['CHROM'] = record['CHROM']
    output['POS'] = record['POS']
    output['REF'] = record['REF']
    output['ALT'] = record['ALT'][0]
    output['VTE'] = record['INFO']['VTE'][0]

def process_vte_sample(record, sample, output):
    output['SAMPLE'] = sample
    sample_map = output['_sample_map']
    output['ROLE'], output['TYPE'] = sample_map[sample]
    output['SCORE'] = record['SAMPLES'][sample].get('VTES', ['NA'])[0]

VTE_CORE_RECORD = [process_vte_record]
VTE_CORE_SAMPLE = [process_vte_sample]
VTE_CORE_HEADER = [process_vte_header]

COSMIC_FIELDS = ['COSMIC_HITTYPE', 'COSMIC_P', 'COSMIC_C', 'COSMIC_COUNT']
def process_cosmic(record, output):
    cos = record['INFO'].get('COSMIC', None)
    if cos is None or cos == ['NONE']:
        ht = p = c = co = 'NA'
    else:
        _cos = []
        for _c in cos:
            ht, p, c, co = _c.split('|')
            _cos.append((ht, p, c, co))
        _cos = sorted(_cos, key = lambda x: int(x[3]))
        ht, p, c, co = _cos[-1]
    output.update(dict(zip(COSMIC_FIELDS, (ht, p, c, co))))

COSMIC_FUNS_RECORD = [process_cosmic]
COSMIC_FUNS_SAMPLE = COSMIC_FUNS_HEADER = []

##
# extract csq format information from a vcf info description
def get_csq_format(csq_descriptor):
    _re_csq = re.search(re_csqformat, csq_descriptor)
    if _re_csq is None:
        raise ValueError("Format search failed")
    return _re_csq.group(1).split('|')


VEP_FIELDS = ['SYMBOL', 'GENE', 'CDNA', 'AA', 'CONSEQUENCE', 'SIFT', 'POLYPHEN', 'CLINSIG', 'DOMAIN', 'MOTIF', 'EXISTINGVAR', 'GMAF']

effect_priority = { # modified from https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl
        'transcript_ablation': 1, # A feature ablation whereby the deleted region includes a transcript feature
        'exon_loss_variant': 1, # A sequence variant whereby an exon is lost from the transcript
        'splice_donor_variant': 2, # A splice variant that changes the 2 base region at the 5' end of an intron
        'splice_acceptor_variant': 2, # A splice variant that changes the 2 base region at the 3' end of an intron
        'stop_gained': 3, # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
        'frameshift_variant': 3, # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
        'stop_lost': 3, # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
        'start_lost': 4, # A codon variant that changes at least one base of the canonical start codon
        'initiator_codon_variant': 4, # A codon variant that changes at least one base of the first codon of a transcript
        'disruptive_inframe_insertion': 5, # An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
        'disruptive_inframe_deletion': 5, # An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
        'inframe_insertion': 5, # An inframe non synonymous variant that inserts bases into the coding sequence
        'inframe_deletion': 5, # An inframe non synonymous variant that deletes bases from the coding sequence
        'missense_variant': 6, # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
        'conservative_missense_variant': 6, # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
        'rare_amino_acid_variant': 6, # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
        'transcript_amplification': 7, # A feature amplification of a region containing a transcript
        'stop_retained_variant': 8, # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
        'synonymous_variant': 8, # A sequence variant where there is no resulting change to the encoded amino acid
        'splice_region_variant': 9, # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
        'incomplete_terminal_codon_variant': 10, # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
        'protein_altering_variant': 11, # A sequence variant which is predicted to change the protein encoded in the coding sequence
        'coding_sequence_variant': 11, # A sequence variant that changes the coding sequence
        'mature_miRNA_variant': 11, # A transcript variant located with the sequence of the mature miRNA
        'exon_variant': 11, # A sequence variant that changes exon sequence
        '5_prime_UTR_variant': 12, # A UTR variant of the 5' UTR
        '5_prime_UTR_premature_start_codon_gain_variant': 12, # snpEff-specific effect, creating a start codon in 5' UTR
        '3_prime_UTR_variant': 12, # A UTR variant of the 3' UTR
        'non_coding_exon_variant': 13, # A sequence variant that changes non-coding exon sequence
        'non_coding_transcript_exon_variant': 13, # snpEff-specific synonym for non_coding_exon_variant
        'non_coding_transcript_variant': 14, # A transcript variant of a non coding RNA gene
        'nc_transcript_variant': 14, # A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
        'intron_variant': 14, # A transcript variant occurring within an intron
        'intragenic_variant': 14, # A variant that occurs within a gene but falls outside of all transcript features. This occurs when alternate transcripts of a gene do not share overlapping sequence
        'INTRAGENIC': 14, # snpEff-specific synonym of intragenic_variant
        'NMD_transcript_variant': 15, # A variant in a transcript that is the target of NMD
        'upstream_gene_variant': 16, # A sequence variant located 5' of a gene
        'downstream_gene_variant': 16, # A sequence variant located 3' of a gene
        'TFBS_ablation': 17, # A feature ablation whereby the deleted region includes a transcription factor binding site
        'TFBS_amplification': 17, # A feature amplification of a region containing a transcription factor binding site
        'TF_binding_site_variant': 17, # A sequence variant located within a transcription factor binding site
        'regulatory_region_ablation': 17, # A feature ablation whereby the deleted region includes a regulatory region
        'regulatory_region_amplification': 17, # A feature amplification of a region containing a regulatory region
        'regulatory_region_variant': 17, # A sequence variant located within a regulatory region
        'regulatory_region':17, # snpEff-specific effect that should really be regulatory_region_variant
        'feature_elongation': 18, # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
        'feature_truncation': 18, # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
        'intergenic_variant': 19, # A sequence variant located in the intergenic region, between genes
        'intergenic_region': 19, # snpEff-specific effect that should really be intergenic_variant
        '': 20
}

biotype_priority = { # modified from https://github.com/mskcc/vcf2maf/blob/master/vcf2maf.pl 
        'protein_coding': 1, # Contains an open reading frame (ORF)
        'LRG_gene': 2, # Gene in a "Locus Reference Genomic" region known to have disease-related sequence variations
        'IG_C_gene': 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_D_gene': 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_J_gene': 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_LV_gene': 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'IG_V_gene': 2, # Immunoglobulin (Ig) variable chain genes imported or annotated according to the IMGT
        'TR_C_gene': 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_D_gene': 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_J_gene': 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'TR_V_gene': 2, # T-cell receptor (TcR) genes imported or annotated according to the IMGT
        'miRNA': 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'snRNA': 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'snoRNA': 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'ribozyme': 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'sRNA': 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'scaRNA': 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'rRNA': 3, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'lincRNA': 3, # Long, intervening noncoding (linc) RNAs, that can be found in evolutionarily conserved, intergenic regions
        'known_ncrna': 4,
        'vaultRNA': 4, # Short non coding RNA genes that form part of the vault ribonucleoprotein complex
        'macro_lncRNA': 4, # unspliced lncRNAs that are several kb in size
        'Mt_tRNA': 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'Mt_rRNA': 4, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'antisense': 5, # Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand
        'sense_intronic': 5, # Long non-coding transcript in introns of a coding gene that does not overlap any exons
        'sense_overlapping': 5, # Long non-coding transcript that contains a coding gene in its intron on the same strand
        '3prime_overlapping_ncrna': 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
        'misc_RNA': 5, # Non-coding RNA predicted using sequences from RFAM and miRBase
        'non_coding': 5, # Transcript which is known from the literature to not be protein coding
        'regulatory_region': 6, # A region of sequence that is involved in the control of a biological process
        'disrupted_domain': 6, # Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain
        'processed_transcript': 6, # Doesn't contain an ORF
        'TEC': 6, # To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies
        'TF_binding_site': 7, # A region of a nucleotide molecule that binds a Transcription Factor or Transcription Factor complex
        'CTCF_binding_site':7, # A transcription factor binding site with consensus sequence CCGCGNGGNGGCAG, bound by CCCTF-binding factor
        'promoter_flanking_region': 7, # A region immediately adjacent to a promoter which may or may not contain transcription factor binding sites
        'enhancer': 7, # A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter
        'promoter': 7, # A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery
        'open_chromatin_region': 7, # A DNA sequence that in the normal state of the chromosome corresponds to an unfolded, un-complexed stretch of double-stranded DNA
        'retained_intron': 7, # Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants
        'nonsense_mediated_decay': 7, # If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD
        'non_stop_decay': 7, # Transcripts that have polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation
        'ambiguous_orf': 7, # Transcript believed to be protein coding, but with more than one possible open reading frame
        'pseudogene': 8, # Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following
        'processed_pseudogene': 8, # Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome
        'polymorphic_pseudogene': 8, # Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated
        'retrotransposed': 8, # Pseudogene owing to a reverse transcribed and re-inserted sequence
        'translated_processed_pseudogene': 8, # Pseudogenes that have mass spec data suggesting that they are also translated
        'translated_unprocessed_pseudogene': 8, # Pseudogenes that have mass spec data suggesting that they are also translated
        'transcribed_processed_pseudogene': 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
        'transcribed_unprocessed_pseudogene': 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
        'transcribed_unitary_pseudogene': 8, #Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
        'unitary_pseudogene': 8, # A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species
        'unprocessed_pseudogene': 8, # Pseudogene that can contain introns since produced by gene duplication
        'Mt_tRNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'tRNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'snoRNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'snRNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'scRNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'rRNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'misc_RNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'miRNA_pseudogene': 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
        'IG_C_pseudogene': 8, # Inactivated immunoglobulin gene
        'IG_D_pseudogene': 8, # Inactivated immunoglobulin gene
        'IG_J_pseudogene': 8, # Inactivated immunoglobulin gene
        'IG_V_pseudogene': 8, # Inactivated immunoglobulin gene
        'TR_J_pseudogene': 8, # Inactivated immunoglobulin gene
        'TR_V_pseudogene': 8, # Inactivated immunoglobulin gene
        'artifact': 9, # Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
        '': 9
}

def parse_vep_header(header, output):
    csq = [h for h in header.get_headers('INFO', 'CSQ')]
    if len(csq) != 1:
        raise ValueError("There were %s CSQ tags, exactly 1 is required", str(len(csq)))
    output['_csqsplit'] = get_csq_format(csq[0]['Description'])

def get_g_pos(a):
    try:
        if a:
            if '-' in a:
                a, b = a.split('-', 1)
            return int(a)
        else:
            return 0
    except ValueError:
        return 0

def cmp_csq(a, b):
    # compare the biotype
    aa = biotype_priority.get(a['BIOTYPE'], 50)
    bb = biotype_priority.get(b['BIOTYPE'], 50)
    if aa != bb:
        return aa - bb
    # compare the effect
    aa = effect_priority.get(a['Consequence'], 9)
    bb = effect_priority.get(b['Consequence'], 9)
    if aa != bb:
        return aa - bb
    aa = a['CLIN_SIG']
    bb = b['CLIN_SIG']
    if aa != bb:
        if aa != '': # means aa must have something
            return -1
        if bb != '': # means bb must have something
            return 1
    try:
        result = get_g_pos(b['cDNA_position']) - get_g_pos(a['cDNA_position']) # longer one gets higher value
        return result
    except:
        print a, b
        raise

def _csq_alleles(record):
    match_alleles = []
    ref = record['REF']
    ref_len = len(ref)
    for a in record['ALT']:
        if ref_len == len(a):
            match_alleles.append(a)
        elif len(a) == 1:
            match_alleles.append('-')
        else:
            match_alleles.append(a[1:])
    return match_alleles

def parse_vep_record(record, output):
    csq_split = output['_csqsplit']
    csq_infos = [dict(zip(csq_split, c.split('|'))) for c in record['INFO'].get('CSQ', [])]
    alts = _csq_alleles(record)    
    csq_infos = sorted([c for c in csq_infos if c['Allele'] in alts], cmp = cmp_csq)
    if len(csq_infos) > 0:
        csq_info = csq_infos[0]
        for k, v in [('MOTIF', 'MOTIF_NAME'), ('DOMAIN', 'DOMAINS'), ('CLINSIG', 'CLIN_SIG'),
                ('SIFT', 'SIFT'), ('POLYPHEN', 'PolyPhen'), ('SYMBOL', 'SYMBOL'), ('GENE', 'Gene'),
                ('CONSEQUENCE', 'Consequence'), ('EXISTINGVAR', 'Existing_variation'),
                ('CDNA', 'HGVSc'), ('AA', 'HGVSp'), ('GMAF', 'GMAF')]:
            output[k] = csq_info[v]
    else:
        for k in VEP_FIELDS:
            output[k] = 'NA'

VEP_FUNS_SAMPLE = []
VEP_FUNS_HEADER = [parse_vep_header]
VEP_FUNS_RECORD = [parse_vep_record]

CAN_FIELDS = ['VALSTATUS', 'CONTEXT_START', 'CONTEXT']

def parse_can_record(record, output):
    output['VALSTATUS'] = ';'.join(record['INFO'].get('DBVS', []))
    output['CONTEXT'] = '|'.join(record['INFO'].get('CONTEXT', []))
    output['CONTEXT_START'] = record['POS'] - output['_context_size']

re_context_size = re.compile(r'''POS - (\d*),''')

def parse_can_header(header, output):
    context_header = [h for h in header.get_headers('INFO', 'CONTEXT')]
    if len(context_header) != 1:
        raise ValueError("More or less than one INFO CONTEXT found: %s" % str(context_header))
    context_size = int(re_context_size.search(context_header[0]['Description']).group(1))
    output['_context_size'] = context_size

CAN_FUNS_SAMPLE = []
CAN_FUNS_HEADER = [parse_can_header]
CAN_FUNS_RECORD = [parse_can_record]

WJ_FIELDS = ['ALLELECOUNT', 'ALLELEQ20COUNT', 'REFCOUNT', 'DP', 'MIDREADCOUNT', 'FORWARDCOUNT', 'REVERSECOUNT', 'CALLERS']

def parse_wj_record(record, output):
    output['_altindex'] = hgsc_vcf.best_alt_index(record)
    output['_refindex'] = hgsc_vcf.ref_index(record)
    oc = [c for o in record['INFO']['OC'] for c in o.split('|')]
    output['CALLERS'] = ';'.join(['/'.join(o.split('/')[-3:]) for o in oc])

def parse_wj_sample(record, sample, output):
    sampleinfo = record['SAMPLES'][sample]
    alt_i = output['_altindex']
    if alt_i < 0:
        output.update({k:'0' for k in ['ALLELECOUNT', 'ALLELEQ20COUNT', 'MIDREADCOUNT', 'FORWARDCOUNT', 'REVERSECOUNT']})
    else:
        for k, v in zip(['ALLELECOUNT', 'ALLELEQ20COUNT', 'MIDREADCOUNT', 'FORWARDCOUNT', 'REVERSECOUNT'] , ['AC', 'AQC', 'MR', 'FC', 'RC']):
            output[k] = sampleinfo[v][alt_i]
    output['DP'] = sampleinfo['DP'][0]
    output['REFCOUNT'] = sampleinfo['AC'][output['_refindex']]

        
WG_FUNS_SAMPLE = [parse_wj_sample]
WJ_FUNS_HEADER = []
WJ_FUNS_RECORD = [parse_wj_record]

def simplify_allele(record):
    record['REF'], alt, record['POS'] = hgsc_vcf._simplify_allele(record['REF'], [record['ALT']], record['POS'])
    record['ALT'] = alt[0]

GATK_FIELDS = ['GT', 'ALLELECOUNT', 'DP']

def parse_gatk_sample(record, sample, output):
    sampleinfo = record['SAMPLES'][sample]
    output['GT'] = sampleinfo['GT'][0]
    output['ALLELECOUNT'] = sampleinfo['AD'][1]
    output['DP'] = sampleinfo['DP'][0]

GATK_FUNS_SAMPLE = [parse_gatk_sample]
GATK_FUNS_HEADER = GATK_FUNS_RECORD = []

def main(args):
    reader = hgsc_vcf.Reader(args.INPUT)
    
    ## 
    # The field names and processing functions must be established based on the command options
    fields = []
    fields += STANDARD_FIELDS
    process_functions_sample = VTE_CORE_SAMPLE
    process_functions_header = VTE_CORE_HEADER
    process_functions_record = VTE_CORE_RECORD
    for _switch, _fields, _sfuns, _hfuns, _rfuns in [
            (args.vep, VEP_FIELDS, VEP_FUNS_SAMPLE, VEP_FUNS_HEADER, VEP_FUNS_RECORD),
            (args.cosmic, COSMIC_FIELDS, COSMIC_FUNS_SAMPLE, COSMIC_FUNS_HEADER, COSMIC_FUNS_RECORD),
            (args.canannot, CAN_FIELDS, CAN_FUNS_SAMPLE, CAN_FUNS_HEADER, CAN_FUNS_RECORD),
            (args.wheeljack, WJ_FIELDS, WG_FUNS_SAMPLE, WJ_FUNS_HEADER, WJ_FUNS_RECORD),
            (args.gatk, GATK_FIELDS, GATK_FUNS_SAMPLE, GATK_FUNS_HEADER, GATK_FUNS_RECORD)
            ]:
        if _switch:
            fields += _fields
            process_functions_sample += _sfuns
            process_functions_header += _hfuns
            process_functions_record += _rfuns
    with bz2.BZ2File(args.OUTPUT, 'w') as fo:
        writer = csv.DictWriter(fo, delimiter = '\t', fieldnames = fields, extrasaction = 'ignore')
        writer.writeheader()
    
        _output_header = {}
        for _f in process_functions_header:
            _f(reader.header, _output_header)
        if args.subject:
            _output_header['SUBJECT'] = args.subject
    
        for record in reader:
            if len(record['REF']) > 50:
                continue
            if 'NONE' in record['INFO']['VTE'][0]:
                continue
            output = {}
            output.update(_output_header)
            for _f in process_functions_record:
                _f(record, output)
            simplify_allele(output)
            for sample in record['SAMPLES']:
                sample_output = {}
                sample_output.update(output)
                for _f in process_functions_sample:
                    _f(record, sample, sample_output)
                for k, v in sample_output.items():
                    if not v:
                        sample_output[k] = 'NA'
                writer.writerow(sample_output)
    

if __name__ == '__main__':
    import argparse
    logger.setLevel(logging.INFO) # adjust the log level to info

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--subject', type = str, help = 'Override the subject id from the VCF, might not be a good idea to do this, but sometimes we have to')
    parser.add_argument('--vep', action = 'store_true', help = 'Activate VEP annotation processing')
    parser.add_argument('--cosmic', action = 'store_true', help = 'Activate COSMIC annotation processing')
    parser.add_argument('--canannot', action = 'store_true', help = 'Activate CANANNOTATION processing')
    parser.add_argument('--wheeljack', action = 'store_true', help = 'Activate Wheeljack processing')
    parser.add_argument('--gatk', action = 'store_true', help = 'Activate GATK processing')
    parser.add_argument('INPUT', type = argparse.FileType('r'), help = 'Input file stream')
    parser.add_argument('OUTPUT', type = str, help = 'Output file stream, will be a bzip file')

    args = parser.parse_args()

    main(args)
