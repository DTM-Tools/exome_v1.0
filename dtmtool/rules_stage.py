import logging
import sys
import traceback
import json
from allele_info import *
from indel_info import *
from multi_info import *
from utils import *
from param_stage import *
from rylanexceptions import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

FILTERS_DEF_JSON = """
    [
        {"name" : "MAPQ","type" : "float","cardinality" : "single"},
        {"name" : "DP","type" : "int","cardinality" : "single"},
        {"name" : "QUAL","type" : "float","cardinality" : "single"},
        {"name" : "AO","type" : "int","cardinality" : "array"},
        {"name" : "AO/DP","type" : "float","cardinality" : "array"}
    ]
    """

INDEL_FILTERS_DEF_JSON = """
    [
        {"name" : "MQM","type" : "float","cardinality" : "array"},
        {"name" : "MQMR","type" : "float","cardinality" : "single"},
        {"name" : "DP","type" : "int","cardinality" : "single"},
        {"name" : "QUAL","type" : "float","cardinality" : "single"},
        {"name" : "AO","type" : "int","cardinality" : "array"},
        {"name" : "AO/DP","type" : "float","cardinality" : "array"},
        {"name" : "QA","type" : "float","cardinality" : "array"},
        {"name" : "QR","type" :  "float","cardinality" : "single"}
    ]
"""

class RylanRules:

    @staticmethod
    def process_rules_stage(pipeline_state):
        rylan_params = pipeline_state['params']
        rylan_rules = RylanRules(rylan_params.get_reference_directory(),
                                 rylan_params.get_rules_file(),
                                 rylan_params.get_rules_indel_file(),
                                 rylan_params.get_rules_multi_file())
        pipeline_state['rules'] = rylan_rules

    def __init__(self, ref_dir, rules_file, rules_indel_file, rules_multi_file):
        self.rules_file = '%s/%s' % (ref_dir, rules_file)
        self.rules_indel_file = '%s/%s' % (ref_dir, rules_indel_file)
        self.rules_multi_file = '%s/%s' % (ref_dir, rules_multi_file)

        self.load_rylan_rules()

    def create_rules_mappings(self):

        self.rules_by_allele = dict()
        self.rules_by_chromo = dict()
        self.spdi_list_hg19 = list()
        self.spdi_list_grch38 = list()
        for a_rule in self.rules_info:
            self.rules_by_allele[a_rule['allele']] = a_rule

            for a_spdi in a_rule['hg19']['spdi']:
                self.spdi_list_hg19.append(a_spdi)

            for a_spdi in a_rule['grch38']['spdi']:
                self.spdi_list_grch38.append(a_spdi)

            if self.rules_by_chromo.get(a_rule['chromosome'], None)== None:
                self.rules_by_chromo[a_rule['chromosome']] = dict()
            self.rules_by_chromo[a_rule['chromosome']][a_rule['allele']] = a_rule

        self.rules_by_allele = sort_dict(self.rules_by_allele)

        self.rules_filters_by_allele = dict()
        for a_filter in self.rules_filters:
            self.rules_filters_by_allele[a_filter['allele']] = a_filter

        self.rules_filters_by_allele = sort_dict(self.rules_filters_by_allele)

        self.indel_by_allele = dict()
        self.indel_by_chromo = dict()
        self.spdi_list_indel_hg19 = list()
        self.spdi_list_indel_grch38 = list()
        for a_indel_rule in self.indel_info:
            self.indel_by_allele[a_indel_rule['allele']] = a_indel_rule

            for a_spdi in a_indel_rule['hg19']['spdi']:
                self.spdi_list_indel_hg19.append(a_spdi)

            for a_spdi in a_indel_rule['grch38']['spdi']:
                self.spdi_list_indel_grch38.append(a_spdi)

            if self.indel_by_chromo.get(a_indel_rule['chromosome'], None) == None:
                self.indel_by_chromo[a_indel_rule['chromosome']] = dict()
            self.indel_by_chromo[a_indel_rule['chromosome']][a_indel_rule['allele']] = a_indel_rule

        self.indel_by_allele = sort_dict(self.indel_by_allele)

        self.indel_filters_by_allele = dict()
        for a_indel_filter in self.indel_filters:
            self.indel_filters_by_allele[a_indel_filter['allele']] = a_indel_filter

        self.indel_filters_by_allele = sort_dict(self.indel_filters_by_allele)

        filters_defs_list = json.loads(FILTERS_DEF_JSON)

        indel_filters_defs_list = json.loads(INDEL_FILTERS_DEF_JSON)

        self.filters_def = dict()
        for a_filter in filters_defs_list:
            self.filters_def[a_filter['name']] = a_filter

        self.indel_filters_def = dict()
        for a_filter in indel_filters_defs_list:
            self.indel_filters_def[a_filter['name']] = a_filter

    def load_rylan_rules(self):

        logger.info('Beggining to read RyLAN rules from files: %s, %s' % (self.rules_file, self.rules_indel_file))

        try:
            self.rules_info, self.rules_filters = read_rylan_allele_info(self.rules_file)
            self.indel_info, self.indel_filters = read_rylan_indel_info(self.rules_indel_file)

            self.create_rules_mappings()

            self.multi_rules = read_multi_info(self.rules_multi_file)

            logger.info('Finished reading RyLAN rules.')

        except Exception as e:
            logger.error('Failed to process rule stage: %s' % str(e))
            traceback.print_exc(file=sys.stdout)
            sys.exit(1)

    def get_spdi_list(self, fasta_ref):
        if (fasta_ref.lower() == 'hg19'):
            return self.spdi_list_hg19
        elif (fasta_ref.lower() == 'grch38'):
            return self.spdi_list_grch38
        else:
            raise new ('Invalid FASTA reference.')

    def get_spdi_list_indel(self, fasta_ref):
        if (fasta_ref.lower() == 'hg19'):
            print(json.dumps(self.spdi_list_indel_hg19))
            return self.spdi_list_indel_hg19
        elif (fasta_ref.lower() == 'grch38'):
            print(json.dumps(self.spdi_list_indel_grch38))
            return self.spdi_list_indel_grch38
        else:
            raise new ('Invalid FASTA reference.')

    def get_multi_rules(self):
        return self.multi_rules

    def get_all_rylan_rules(self):
        return self.rules_by_allele

    def get_rylan_rule(self, allele):
        return self.rules_by_allele.get(allele, None)

    def get_all_rylan_rules_by_chromo(self):
        return self.rules_by_chromo

    def get_rylan_rules_by_chromo(self, chromo_id):
        return self.rules_by_chromo.get(chromo_id, None)

    def get_all_rylan_filters(self):
        return self.rules_filters_by_allele

    def get_rylan_filter(self, allele):
        return self.rules_filters_by_allele.get(allele, None)

    def get_all_rylan_indel_rules(self):
        return self.indel_by_allele

    def get_rylan_indel_rule(self, allele):
        return self.indel_by_allele.get(allele, None)

    def get_all_rylan_indel_rules_by_chromo(self):
        return self.indel_by_chromo

    def get_rylan_indel_rules_by_chromo(self, chromo_id):
        return self.indel_by_chromo.get(chromo_id, None)

    def get_all_rylan_indel_filters(self):
        return self.indel_filters_by_allele

    def get_rylan_indel_filter(self, allele):
        return self.indel_filters_by_allele.get(allele, None)

    def get_filter_definitions(self):
        return self.filters_def

    def get_indel_filter_definitions(self):
        return self.indel_filters_def
