import os
import sys
import multiprocessing
import subprocess
import datetime
import logging
import json

from utils import *
from rules_stage import *
from param_stage import *
from vcf_rylan_stage import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class RylanInDelVcfData:

    @staticmethod
    def process_vcf_indel_stage(pipeline_state):
        vcf_indel_container = RylanInDelVcfData(pipeline_state['params'], pipeline_state['rules'], pipeline_state['vcfInfo'])
        vcf_indel_container.generate_vcfs_indel()
        vcf_indel_container.evaluate_vcf_files_indel()


    def __init__(self, rylan_params, rylan_rules, rylan_info):
        self.rylan_params = rylan_params
        self.rylan_rules = rylan_rules
        self.rylan_info = rylan_info

    def generate_vcfs_indel(self):
        indel_classifiers_by_chromo = self.rylan_rules.get_all_rylan_indel_rules_by_chromo()
        ref_genome = self.rylan_params.get_reference_build()
        ref_dir = self.rylan_params.get_reference_directory()
        input_dir = self.rylan_params.get_input_directory()
        patient_file = self.rylan_params.get_patient_file()

        indel_ranges_by_chromo = RylanInDelVcfData.get_chromo_positions_info_indel(indel_classifiers_by_chromo, ref_genome)

        pool = multiprocessing.Pool(self.rylan_params.get_num_cores())

        bayes_args=[]
        for chromo in indel_classifiers_by_chromo.keys():
            an_arg = (self.rylan_info.get_run_id(), ref_dir, input_dir, self.rylan_info.get_temp_dir(), ref_genome,
                    patient_file, chromo, indel_ranges_by_chromo[chromo])
            bayes_args.append(an_arg)

        logger.info('Run ID %s:Starting freebayes parallel execution ...' % self.rylan_info.get_run_id())
        temp_freebayes_res = pool.map(RylanVcfData.execute_freebayes_mp, bayes_args)
        logger.info('Run ID %s:Completed FreeBayes execution ...' % self.rylan_info.get_run_id())

        logger.info('Run ID %s: Computing MAPQ scores ...' % self.rylan_info.get_run_id())
        indel_ranges_by_chromo_no_flatten = RylanInDelVcfData.get_chromo_positions_info_indel_no_flatten(indel_classifiers_by_chromo, ref_genome)

        for chromo in indel_ranges_by_chromo.keys():
            mapq_scores = RylanVcfData.determine_mapq(input_dir, patient_file, chromo,
                indel_ranges_by_chromo_no_flatten[chromo])
            self.rylan_info.get_mapq_results().update(mapq_scores)

        return


    @staticmethod
    def read_vcf_file_indel(run_id, vcf_file_name):
        logger.info('Run ID %s: Reading VCF file: %s' % (run_id, vcf_file_name))
        fdat = open(vcf_file_name, 'r')
        vcf_data = vcf.Reader(fdat)
        records = []

        try:
            for record in vcf_data:
                records.append(record)

            fdat.close()
            return records
        except:
            fdat.close()
            raise RylanProcessException('unexpected: error reading VCF file... Abort.')

    @staticmethod
    def apply_allele_filter_indel(in_result, rules, filter_defs):

        in_result['filter_criteria'] = []
        rule_set = []
        if in_result['determination'] == 'classified':
            if in_result['variant_count'] == 0:

                rule_set = rules['reference']
            elif in_result['zygosity'] == 1:
                rule_set = rules['variant']
            else:
                rule_set = rules['homozygous_variant']

            for rule in rule_set:
                filter_info = filter_defs[rule['filterId']]
                if filter_info['cardinality'] == 'single':
                    filter_eval = RylanVcfData.evaluate_filter(rule['criteria'], in_result[rule['filterId']])
                    if filter_eval == False:
                        if rule['action'] == 'flag':

                            in_result['determination'] = 'filtered'
                            in_result['filter_criteria'].append('%s: %s'%(rule['filterId'], rule['criteria']))
                        elif rule['action'] == 'warn':
                            in_result['determination'] = 'unverified'
                            in_result['filter_criteria'].append('%s: %s'%(rule['filterId'], rule['criteria']))
                elif filter_info['cardinality'] == 'array':
                    for param_value in in_result[rule['filterId']]:
                        filter_eval = RylanVcfData.evaluate_filter(rule['criteria'], param_value)
                        if filter_eval == False:
                            if rule['action'] == 'flag':
                            
                                in_result['determination'] = 'filtered'
                                in_result['filter_criteria'].append('%s: %s'%(rule['filterId'], rule['criteria']))
                            elif rule['action'] == 'warn':
                                in_result['determination'] = 'unverified'
                                in_result['filter_criteria'].append('%s: %s'%(rule['filterId'], rule['criteria']))

                else:
                    raise Exception("Error in filter definition")
            return in_result
        else:
            return in_result


    def evaluate_vcf_files_indel(self):
        ref_genome = self.rylan_params.get_reference_build()
        filter_rules = self.rylan_rules.get_all_rylan_indel_filters()
        filter_defs = self.rylan_rules.get_indel_filter_definitions()
        tmp_dir = self.rylan_info.get_temp_dir()
        run_id = self.rylan_info.get_run_id()
        filtered_results = self.rylan_info.get_filtered_results()
        mapq_results = self.rylan_info.get_mapq_results()

        for chromo in CHROMO_LIST:
            chromo_classifiers = RylanInDelVcfData.get_enabled_indel_rules_by_chromosome(self.rylan_rules.get_all_rylan_indel_rules(), chromo)
            for allele in chromo_classifiers.keys():
                print(">>> Processing allele: %s"%(allele))
                rule_data = chromo_classifiers[allele]
                vcf_file_name = "%s%s.%s.%s.vcf"%(tmp_dir, run_id, chromo, allele.replace("/","_"))


                vcf_records = RylanInDelVcfData.read_vcf_file_indel(run_id, vcf_file_name)
                vcf_data = None
                if len(vcf_records) > 1:
                    logger.warn('More than one VCF record read: %d'%len(vcf_records))
                    filtered_results[allele] = RylanVcfData.mark_allele_negative(allele, rule_data)
                elif len(vcf_records) == 0:
                    logger.warn('Run ID %s: Allele not found: %s' % (run_id, allele))
                    filtered_results[allele] = RylanVcfData.mark_allele_not_found(allele, rule_data)

                elif vcf_records == None:
                    logger.warn('Run ID %s: Allele not found: %s' % (run_id, allele))
                    raise RylanProcessException('unexpected: more than one record per VCF file... Abort.')
                else:

                    vcf_data = vcf_records[0]
                    if (vcf_data != None):
                        unfiltered_result = RylanInDelVcfData.evaluate_vcf_rules_indel(vcf_data, rule_data, ref_genome, mapq_results)
                        filtered_result = RylanInDelVcfData.apply_allele_filter_indel(unfiltered_result, filter_rules[allele], filter_defs)
                        filtered_results[allele] = filtered_result
                    else:
                        logger.warn('Run ID %s: Allele not found: %s' % (run_id, allele))
                        filtered_results[allele] = RylanVcfData.mark_allele_not_found(allele, rule_data)


        return

    @staticmethod
    def evaluate_vcf_rules_indel(record, rule_data, fasta_ref, mapq_results):
        rule_variants = rule_data[fasta_ref]
        print(rule_variants)
        print(record)
        info = {}

        info['allele'] = rule_data['allele']
        info['description'] = rule_data['description']
        info['position_start'] = int(rule_variants['position_start'])
        info['position_end'] = int(rule_variants['position_end'])
        info['cell'] = rule_data['cell']
        info['alleleType'] = rule_data['function']
        info['genSubtype'] = rule_data['genSubtype']

        if len(record.ALT) == 1 and str(record.ALT[0]) == '<*>':

            info['determination'] = 'classified'
            info['classification'] = rule_variants['reference']['class']
            info['zygosity'] = 0
            info['found'] = [str(record.REF)]
            info['effect'] = rule_data['refSubtype']
            info['variant_count'] = 0
            info['DP'] = record.INFO['DP']
            info['QR'] = record.samples[0]['QR']
            info['QA'] = [record.samples[0]['QA']]
            info['MQM'] = [100000]
            info['MQMR'] = 100000
            info['MAPQ'] = mapq_results[rule_data['allele']]
            if (record.is_indel):
                if (record.is_deletion):
                    info['type'] = 'deletion'
                else:
                    info['type'] = 'insertion'
            else:
                info['type'] = 'normal'

            return info
        elif len(record.ALT) == 1 and str(record.ALT[0]) != '<*>':

            info['found'] = [str(record.ALT[0])]
            alternate = rule_variants['alternate']
            if alternate['nucleotides'] == str(record.ALT[0]):
                info['determination'] = 'classified'
                info['variant_count'] = 1
                info['zygosity'] = record.samples[0].gt_type
                if record.samples[0].gt_type == 1:
                    info['classification'] = '%s/%s'%(rule_variants['reference']['class'],alternate['class'])
                    info['effect'] = rule_data['refSubtype'],rule_data['altSubtype']
                elif record.samples[0].gt_type == 0:
                    info['classification'] = rule_variants['reference']['class']
                    info['effect'] = rule_data['refSubtype']
                else:
                    info['classification'] = alternate['class']
                    info['effect'] = rule_data['altSubtype']
                info['QUAL'] = record.QUAL
                info['DP'] = record.INFO['DP']
                info['AO'] = record.INFO['AO']
                info['AO/DP'] = [info['AO'][0] / info['DP']]
                info['QR'] = record.INFO['QR']
                info['QA'] = record.INFO['QA']
                info['MQM'] = record.INFO['MQM']
                info['MQMR'] = record.INFO['MQMR']
                info['MAPQ'] = mapq_results[rule_data['allele']]
                if (record.is_indel):
                    if (record.is_deletion):
                        info['type'] = 'deletion'
                    else:
                        info['type'] = 'insertion'
                else:
                    info['type'] = 'normal'
                return info
            else:

                info['determination'] = 'abnormality'
                return info
        elif len(record.ALT) == 2:

            info['found'] = [str(record.ALT[0]), str(record.ALT[1])]

            info['determination'] = 'inconsistent'
            info['expected'] = rule_variants['reference']['nucleotides']
            return info
        else:

            info['determination'] = 'inconsistent'
            info['expected'] = rule_variants['reference']['nucleotides']
            info['found'] = [str(record.ALT[0]),str(record.ALT[1]),str(record.ALT[2])]
            return info

    @staticmethod
    def get_chromo_positions_info_indel_no_flatten(chromo_classifiers, ref_genome):
        ranges_by_chromo = dict()
        for chromo in chromo_classifiers.keys():
            positions_list = RylanInDelVcfData.get_position_list_indel_no_flatten(chromo_classifiers[chromo], ref_genome)
            ranges_by_chromo[chromo] = positions_list
        return ranges_by_chromo

    @staticmethod
    def get_chromo_positions_info_indel(classifier_rules, ref_genome):
        ranges_by_chromo = {}
        for chromo in classifier_rules.keys():
            positions_list = RylanInDelVcfData.get_position_list_indel(classifier_rules[chromo], ref_genome)
            ranges_by_chromo[chromo] = positions_list
        return ranges_by_chromo

    @staticmethod
    def get_position_list_indel(classifier_rules, ref_genome):
        positions_found = []
        for allele in classifier_rules.keys():
            pos_value_st = classifier_rules[allele][ref_genome]['position_start']
            pos_value_end = classifier_rules[allele][ref_genome]['position_end']
            pos_found = (allele.replace('/','_'),pos_value_st, pos_value_end)
            positions_found.append(pos_found)

        return positions_found

    @staticmethod
    def get_position_list_indel_no_flatten(classifier_rules, ref_genome):
        positions_found = []
        for allele in classifier_rules.keys():
            pos_value_st = classifier_rules[allele][ref_genome]['position_start']
            pos_value_end = classifier_rules[allele][ref_genome]['position_end']
            pos_found = (allele,pos_value_st, pos_value_end)
            positions_found.append(pos_found)

        return positions_found

    @staticmethod
    def get_enabled_indel_rules_by_chromosome(classifier_collection, chromo_id):
        classifiers = dict()
        for an_allele in classifier_collection.keys():
            a_rule = classifier_collection[an_allele]

            if (a_rule['enabled'] == True) and (a_rule['chromosome'] == chromo_id):
                classifiers[a_rule['allele']] = a_rule
        return classifiers
