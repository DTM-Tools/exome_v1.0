import sys
import os
import logging
import json

from rylanexceptions import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

COL_ALLELE = 0
COL_DESC = COL_ALLELE + 1
COL_LIST = COL_DESC + 1
COL_REF = COL_LIST + 1
COL_ALT = COL_REF + 1

def read_multi_info(filein):

    try:
        logger.info('Reading Multi definition file...')
        file_var = open(filein,"r")

        lines = list()

        output_rules = dict()

        line_count = 0
        for line in file_var:
            if line_count != 0:
                line = line.replace('\n','').replace('\r','')
                file_vals = line.split(',')
                multi_name = file_vals[COL_ALLELE]
                alleles_desc = file_vals[COL_DESC]
                allele_list_str = file_vals[COL_LIST]
                allele_list = allele_list_str.split('|')
                reference_classif = file_vals[COL_REF]
                alternate_classif = file_vals[COL_ALT]
                rule_summary = dict()
                rule_summary['alleles'] = allele_list
                rule_summary['description'] = alleles_desc
                rule_summary['refclass'] = reference_classif
                rule_summary['altclass'] = alternate_classif
                output_rules[multi_name] = rule_summary
            line_count = line_count + 1

        file_var.close()
        return output_rules

    except Exception as e:
        logger.error(e)
        raise RylanProcessException('Failed compound allele data.')
