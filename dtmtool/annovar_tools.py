import sys
import os
import logging
import json

from rylanexceptions import *
from utils import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

COL_LOCATION = 0
COL_GENE = COL_LOCATION + 1
COL_CHROM = COL_GENE + 1
COL_START = COL_CHROM + 1
COL_END = COL_START + 1
COL_REF = COL_END + 1
COL_ALT = COL_REF + 1
COL_ZYGO = COL_ALT + 1

def read_annovar_info(filein, filein2):

    try:
        logger.info('Reading Annovar variant var_type definition files for SNVs...')
        file_annovar1 = open(filein,"r")
        file_annovar2 = open(filein2,"r")

        findings_annovarSNV = dict()
        count = 1


        annovar_2 = dict()
        for line in file_annovar2:
            line = line.replace('\n','')
            file_vals = line.split('\t')
            line_number = file_vals[0]
            annovar_2[line_number] = file_vals

        file_annovar2.close()

        for line in file_annovar1:
            line = line.replace('\n','')
            file_vals = line.split('\t')
            location = file_vals[COL_LOCATION]
            gene = file_vals[COL_GENE]
            chrom = file_vals[COL_CHROM]
            start = file_vals[COL_START]
            end = file_vals [COL_END]
            reference = file_vals[COL_REF]
            alternate = file_vals [COL_ALT]
            annovar_zygosity = file_vals [COL_ZYGO]
            if annovar_zygosity == 'hom':
                zygo_code = 2
            else:
                zygo_code = 1
            findings_summary = dict()
            findings_summary['location'] = location
            findings_summary['gene'] = gene
            findings_summary['chromosome'] = chrom
            findings_summary['start'] = start
            findings_summary['end'] = end
            findings_summary['reference'] = reference
            findings_summary['alternate'] = alternate
            findings_summary['zygosity'] = zygo_code
            line_str = str('line%d')%(count)
            identifier = calculate_spdi_code(chrom,start,reference,alternate)
            findings_summary['identifier'] = identifier

            if findings_summary['location'] == 'exonic':
                line_info = annovar_2.get(line_str, None)
                if line_info != None:
                    var_type = line_info[1]
                    effect_list_str = line_info[2]
                    effect_list = effect_list_str.split(',')[:-1]
                    findings_summary['varType'] = var_type
                    findings_summary ['codingEffect'] = effect_list

            if findings_summary['location'] == 'UTR3' or findings_summary['location'] == 'UTR5':
                str_gene_list_start = findings_summary['gene'].find('(')
                str_list_to_process = findings_summary['gene'][str_gene_list_start + 1: -1]
                # TODO: add UTR effect list
                list_utr = str_list_to_process.split(',')
                findings_summary['UTREffect'] = list_utr
                findings_summary['gene'] = findings_summary['gene'][0:str_gene_list_start - 1]

            findings_annovarSNV[identifier] = findings_summary

            count = count + 1

        file_annovar1.close()
        return findings_annovarSNV

    except Exception as e:
        logger.error(e)
        raise RylanProcessException('Failed annovar definition files')
