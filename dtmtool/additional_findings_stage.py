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
from annovar_tools import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class RylanAdditionalFindings:

    @staticmethod
    def process_additional_findings_stage(pipeline_state):
        additional_findings_container = RylanAdditionalFindings(pipeline_state['params'], pipeline_state['rules'],pipeline_state['vcfInfo'])
        additional_findings_container.process_additional_findings()
        additional_findings_container.filter_findings()
        pipeline_state['additionalFindings'] = additional_findings_container

    def compute_annovar_file_names(self):
        base_name = self.rylan_params.get_patient_file()[:-10]
        input_dir = self.rylan_params.get_input_directory()
        file2_name = '%s/%s_Annovar.exonic_variant_function'%(input_dir, base_name)
        file1_name = '%s/%s_Annovar.variant_function'%(input_dir, base_name)
        return (file1_name, file2_name)

    def __init__(self, rylan_params, rylan_rules, rylan_info):
        self.rylan_params = rylan_params
        self.rylan_rules = rylan_rules
        self.rylan_info = rylan_info

    def get_annovar_findings(self):
        return self.annovar_findings

    def get_annovar_filtered_findings(self):
        return self.annovar_filtered_findings

    def filter_findings(self):
        self.annovar_filtered_findings = dict()
        spdi_codes = self.rylan_rules.get_spdi_list(self.rylan_params.get_reference_build())
        print(spdi_codes)
        spdi_codes_indel = self.rylan_rules.get_spdi_list_indel(self.rylan_params.get_reference_build())
        print(spdi_codes_indel)

        for a_finding_key in self.annovar_findings.keys():
            if a_finding_key not in spdi_codes and a_finding_key not in spdi_codes_indel:
                self.annovar_filtered_findings[a_finding_key] = self.annovar_findings[a_finding_key]

    def process_additional_findings(self):

        self.rules_by_chromo_and_pos = dict()

        rylan_rules_by_chromo = self.rylan_rules.get_all_rylan_rules_by_chromo()

        rylan_rules= self.rylan_rules.get_all_rylan_rules()

        rylan_rules_indel = self.rylan_rules.get_all_rylan_indel_rules()


        rules_SNV_spdi = list()
        annovar_file_tup = self.compute_annovar_file_names()
        self.annovar_findings = read_annovar_info(annovar_file_tup[0], annovar_file_tup[1])
