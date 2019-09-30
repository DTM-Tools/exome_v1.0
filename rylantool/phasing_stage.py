import os
import sys
import multiprocessing
import subprocess
import datetime
import logging
import json

import vcf

from utils import *
from rules_stage import *
from param_stage import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class RylanPhasing:

    @staticmethod
    def process_phasing_stage(pipeline_state):
        phasing_container = RylanPhasing(pipeline_state['params'], pipeline_state['rules'],pipeline_state['vcfInfo'])
        phasing_container.process_phasing()

    def __init__(self, rylan_params, rylan_rules, rylan_info):
        self.rylan_params = rylan_params
        self.rylan_rules = rylan_rules
        self.rylan_info = rylan_info

    def determine_zygosity_classification(self, allele_list):
        all_twos = True
        ones_count = 0;
        at_least_one_zero = False
        filtered_count = 0
        vcf_data = self.rylan_info.get_filtered_results()
        zygosity_data = dict()

        for allele_id in allele_list:
            allele_info = vcf_data.get(allele_id, None)
            if allele_info != None:
                if allele_info['determination'] != 'classified' and allele_info['determination'] != 'filtered':
                    return (-1, 0, list())

                if allele_info['determination'] == 'filtered':
                    filtered_count = filtered_count + 1
                zygosity = allele_info['zygosity']
                zygosity_data[allele_id] = zygosity
                if zygosity != 2:
                    all_twos = False
                if zygosity == 0:
                    at_least_one_zero = True
                if zygosity == 1:
                    ones_count = ones_count + 1
            else:
                return (-1, 0, list())

        if all_twos == True:
            return (2, 1, zygosity_data, filtered_count) # zygosity = 2, certainty = true
        elif at_least_one_zero == True:
            return (0, 1, zygosity_data, filtered_count)
        elif ones_count == 1:
            return (1, 1, zygosity_data, filtered_count)
        else:
            return (1, 0, zygosity_data, filtered_count)

    def process_phasing(self):
        phasing_rules = self.rylan_rules.get_multi_rules()

        for multi_name in phasing_rules.keys():
            findings_det = self.determine_zygosity_classification(phasing_rules[multi_name]['alleles'])
            logger.info('Multi: %s, determination: %d, certainty:%d'%(multi_name, findings_det[0], findings_det[1]))

            findings_out = dict()
            findings_out['alleleList'] = findings_det[2]
            findings_out['allele'] = multi_name
            findings_out['description'] = phasing_rules[multi_name]['description']
            if (findings_det[0] != -1):
                findings_out['zygosity'] = findings_det[0]
                findings_out['certainty'] = findings_det[1]
                if (findings_det[3] == 0):
                    findings_out['determination'] = 'classified'
                else:
                    findings_out['determination'] = 'filtered'
                    findings_out['filteredCount'] = findings_det[3]

                if findings_det[0] == 0:
                    classif_str = phasing_rules[multi_name]['refclass']
                elif findings_det[0] == 2:
                    classif_str = phasing_rules[multi_name]['altclass']
                elif findings_det[0] == 1:
                    if findings_det[1] == 1:
                        classif_str = '%s|%s'%(phasing_rules[multi_name]['refclass'], phasing_rules[multi_name]['altclass'])
                    else:
                        classif_str = '%s|%s-possible'%(phasing_rules[multi_name]['refclass'], phasing_rules[multi_name]['altclass'])

                findings_out['classification'] = classif_str
            else:
                findings_out['determination'] = 'not classfied'


            rylan_findings = self.rylan_info.get_filtered_results()
            rylan_findings[multi_name] = findings_out
