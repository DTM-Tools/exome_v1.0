
import os
import random, string

from rylanexceptions import *

def sort_dict(in_dict):

    new_dict = dict()

    for key in sorted(in_dict.keys()):
        new_dict[key] = in_dict[key]

    return new_dict

def randomword(length):
   return ''.join(random.sample(string.ascii_lowercase + string.digits, length))

def getFastaFileName(reference):

    if reference.lower() == "hg19":
        return 'human_g1k_v37.fasta'
    elif reference.lower() == "grch38":
        return 'GRCh38.fasta'
    else:
        raise RylanParamException('Invalid reference type: %s' % reference)

def calculate_spdi_code(chrom, start, reference, alternate):
    spdi_code = 'chrom%s:%s:%s:%s'%(chrom,start,reference,alternate)
    return spdi_code.upper()

def calculate_spdi_code_ins(chrom, end, reference, alternate):
    ref_corr = "-"
    alt_corr = alternate[1:-1]
    spdi_code = 'chrom%s:%s:%s:%s'%(chrom,end,ref_corr,alt_corr)
    return spdi_code.upper()

def calculate_spdi_code_del(chrom, start, reference, alternate):
    start_corr = int(start) + 2
    start_corr_str = str(start_corr)
    ref_corr = reference[1:-1]
    alt_corr = '-'
    spdi_code = 'chrom%s:%s:%s:%s'%(chrom,start_corr_str,ref_corr,alt_corr)
    return spdi_code.upper()
