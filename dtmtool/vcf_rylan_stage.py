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


CHROMO_LIST = [ "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" ]


class RylanVcfData:

    @staticmethod
    def process_vcf_stage(pipeline_state):
        vcf_container = RylanVcfData(pipeline_state['params'], pipeline_state['rules'])
        vcf_container.generate_vcfs()
        vcf_container.evaluate_vcf_files()
        pipeline_state['vcfInfo'] = vcf_container

    def __init__(self, rylan_params, rylan_rules):
        self.run_id = randomword(16)
        self.rylan_params = rylan_params
        self.rylan_rules = rylan_rules
        self.mapq_results = dict()
        self.filtered_results = dict()

        try:
            output_dir = rylan_params.get_output_directory()
            self.temp_dir = '%s/tmp/'%(output_dir)
            if (not os.path.exists(self.temp_dir)):
                logger.info('Creating temporary file directory: %s' % self.temp_dir)
                os.makedirs(self.temp_dir)
        except Exception as e:
            logger.error(str(e))
            raise RylanProcessException('Failed to create temporary file directory.')

    def get_run_id(self):
        return self.run_id

    def get_temp_dir(self):
        return self.temp_dir

    def get_mapq_results(self):
        return self.mapq_results

    def get_filtered_results(self):
        return self.filtered_results

    def evaluate_vcf_files(self):
        ref_genome = self.rylan_params.get_reference_build()
        filter_rules = self.rylan_rules.get_all_rylan_filters()
        filter_defs = self.rylan_rules.get_filter_definitions()


        for chromo in CHROMO_LIST:
            chromo_classifiers = RylanVcfData.get_enabled_rules_by_chromosome(self.rylan_rules.get_all_rylan_rules(), chromo)
            for allele in chromo_classifiers.keys():
                logger.info('Run ID %s: Processing allele: %s' % (self.run_id, allele))

                rule_data = chromo_classifiers[allele]
                vcf_file_name = '%s%s.%s.%s.vcf' % (self.temp_dir, self.run_id, chromo, allele.replace("/","_"))
                vcf_data = RylanVcfData.read_vcf_file(self.run_id, vcf_file_name)
                if (vcf_data != None):
                    unfiltered_result = RylanVcfData.evaluate_vcf_rules(vcf_data, rule_data, ref_genome, self.mapq_results)
                    filtered_result = RylanVcfData.apply_allele_filter(unfiltered_result, filter_rules[allele], filter_defs)
                    self.filtered_results[allele] = filtered_result
                else:
                    logger.warn('Run ID %s: Allele not found: %s' % (self.run_id, allele))
                    self.filtered_results[allele] = RylanVcfData.mark_allele_not_found(allele, rule_data)

        return

    def generate_vcfs(self):
        classifiers_by_chromo = self.rylan_rules.get_all_rylan_rules_by_chromo()
        ref_genome = self.rylan_params.get_reference_build()
        ref_dir = self.rylan_params.get_reference_directory()
        input_dir = self.rylan_params.get_input_directory()
        patient_file = self.rylan_params.get_patient_file()

        ranges_by_chromo = RylanVcfData.get_chromo_positions_info(classifiers_by_chromo, ref_genome)

        logger.info('Multiprocessing with cores: %d'%self.rylan_params.get_num_cores())
        pool = multiprocessing.Pool(int(self.rylan_params.get_num_cores()))

        bayes_args=[]
        for chromo in classifiers_by_chromo.keys():
            an_arg = (self.run_id, ref_dir, input_dir, self.temp_dir, ref_genome,
                    patient_file, chromo, ranges_by_chromo[chromo])
            bayes_args.append(an_arg)

        logger.info('Run ID %s:Starting freebayes parallel execution ...' % self.run_id)
        temp_freebayes_res = pool.map(RylanVcfData.execute_freebayes_mp, bayes_args)
        logger.info('Run ID %s:Completed FreeBayes execution ...' % self.run_id)

        logger.info('Run ID %s: Computing MAPQ scores ...' % self.run_id)
        ranges_by_chromo_no_flatten = RylanVcfData.get_chromo_positions_info_no_flatten(classifiers_by_chromo, ref_genome)

        for chromo in ranges_by_chromo.keys():

            mapq_scores = RylanVcfData.determine_mapq(input_dir, patient_file, chromo,
                ranges_by_chromo_no_flatten[chromo])
            self.mapq_results.update(mapq_scores)

        return

    def cleanup_vcf_files(self):
        files = os.listdir(self.temp_dir)
        for file in files:
            if (file.startswith(self.run_id)):
                os.remove(os.path.join(self.temp_dir,file))

    @staticmethod
    def evaluate_vcf_rules(record, rule_data, fasta_ref, mapq_results):
        rule_variants = rule_data[fasta_ref]
        print(rule_variants)
        info = {}

        info['allele'] = rule_data['allele']
        info['description'] = rule_data['description']
        info['position'] = int(rule_variants['position'])
        info['cell'] = rule_data['cell']
        info['alleleType'] = rule_data['type']
        info['genSubtype'] = rule_data['genSubtype']


        if len(record.ALT) == 1 and str(record.ALT[0]) == '<*>':

            if str(record.REF) != rule_variants['reference']['nucleotide']:
                info['determination'] = 'inconsistent'
                info['expected'] = rule_variants['reference']['nucleotide']
                info['found'] = [str(record.REF)]
                info['variant_count'] = 0
                return info
            else:
                info['determination'] = 'classified'
                info['classification'] = rule_variants['reference']['classification']
                info['zygosity'] = 0
                info['found'] = [str(record.REF)]
                info['effect'] = rule_data[str(record.REF).lower() + 'Subtype']
                info['variant_count'] = 0
                info['DP'] = record.INFO['DP']
                info['QR'] = record.samples[0]['QR']
                info['MAPQ'] = mapq_results[rule_data['allele']]
                return info
        elif len(record.ALT) == 1 and str(record.ALT[0]) != '<*>':

            info['found'] = [str(record.ALT[0])]
            for variant in rule_variants['variants']:
                if variant['nucleotide'] == str(record.ALT[0]):
                    info['determination'] = 'classified'
                    info['variant_count'] = 1
                    info['zygosity'] = record.samples[0].gt_type
                    if record.samples[0].gt_type == 1:
                        info['classification'] = '%s/%s'%(rule_variants['reference']['classification'],variant['classification'])
                        info['effect'] = [rule_data[str(record.REF).lower() + 'Subtype'], rule_data[str(record.ALT[0]).lower() + 'Subtype']]
                    else:
                        info['classification'] = variant['classification']
                        info['effect'] = rule_data[str(record.ALT[0]).lower() + 'Subtype']
                    info['QUAL'] = record.QUAL
                    info['DP'] = record.INFO['DP']
                    info['AO'] = record.INFO['AO']
                    info['AO/DP'] = [info['AO'][0] / info['DP']]
                    info['MAPQ'] = mapq_results[rule_data['allele']]
                    return info

            info['determination'] = 'abnormality'
            return info
        elif len(record.ALT) == 2:

            info['found'] = [str(record.ALT[0]), str(record.ALT[1])]
            try:
                info['effect'] = [rule_data[str(record.ALT[0]).lower() + 'Subtype'], rule_data[str(record.ALT[1]).lower() + 'Subtype']]
            except:
                pass

            matches_found = 0
            for nuc in record.ALT:
                for variant in rule_variants['variants']:
                    if variant['nucleotide'] == str(nuc):
                        matches_found = matches_found + 1

            if matches_found == 2:
                info['determination'] = 'classified'
                info['classification'] = '%s/%s'%(rule_variants['variants'][0], rule_variants['variants'][1])
                info['variant_count'] = 2
                info['zygosity'] = 1
                info['QUAL'] = record.QUAL
                info['DP'] = record.INFO['DP']
                info['AO'] = record.INFO['AO']
                info['AO/DP'] = [info['AO'][0] / info['DP'], info['AO'][1] / info['DP']]
                info['MAPQ'] = mapq_results[rule_data['allele']]
                return info
            else:

                info['determination'] = 'abnormality'
                return info
        else:

            info['determination'] = 'triples'
            info['expected'] = rule_variants['reference']['nucleotide']
            info['found'] = [str(record.ALT[0]),str(record.ALT[1]),str(record.ALT[2])]
            try:
                info['effect'] = [rule_data[str(record.ALT[0]).lower() + 'Subtype'], rule_data[str(record.ALT[1]).lower() + 'Subtype'], rule_data[str(record.ALT[2]).lower() + 'Subtype']]
            except:
                pass
            return info


    @staticmethod
    def evaluate_filter(expr, x):
        try:
            return eval(expr)
        except Exception:
            raise RylanProcessException('Error evaluating filter: %s' % (expr))


    @staticmethod
    def apply_allele_filter(in_result, rules, filter_defs):


        in_result['filter_criteria'] = []
        rule_set = []

        if in_result['determination'] == 'classified':
            if in_result['variant_count'] == 0:

                rule_set = rules['reference']
            else:

                rule_set = rules['variant']

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

    @staticmethod
    def mark_allele_not_found(allele, rule_data):
        info = dict()
        info['allele'] = rule_data['allele']
        info['description'] = rule_data['description']
        info['determination'] = 'MISSING'
        return info

    @staticmethod
    def mark_allele_negative(allele, rule_data):
        info = dict()
        info['allele'] = rule_data['allele']
        info['description'] = rule_data['description']
        info['determination'] = 'Negative probably other variants'
        return info

    @staticmethod
    def read_vcf_file(run_id, vcf_file_name):
        logger.info('Run ID %s: Reading VCF file: %s' % (run_id, vcf_file_name))
        fdat = open(vcf_file_name, 'r')
        vcf_data = vcf.Reader(fdat)
        records = []

        try:
            for record in vcf_data:

                records.append(record)

            if len(records) > 1:
                raise RylanProcessException('unexpected: more than one record per VCF file... Abort.')

            fdat.close()
            return records[0]
        except:

            fdat.close()
            return None

    @staticmethod
    def get_enabled_rules_by_chromosome(classifier_collection, chromo_id):
        classifiers = dict()
        for an_allele in classifier_collection.keys():
            a_rule = classifier_collection[an_allele]

            if (a_rule['enabled'] == True) and (a_rule['chromosome'] == chromo_id):
                classifiers[a_rule['allele']] = a_rule
        return classifiers


    @staticmethod
    def determine_mapq(input_dir, patient_file, chromo, pos_ranges):
        rylan_dir = os.environ.get('RYLAN_DIR', '.')
        file_name = '%s/%s'%(input_dir, patient_file)

        mapq_res = dict()
        for one_range in pos_ranges:
            range_arg = '%s:%s-%s'%(chromo, one_range[1], one_range[2])
            exec_args = ['sh', '%s/dtmtool/mapqcalc.sh' % rylan_dir, file_name, range_arg]
            logger.info('executing: %s'%(exec_args))

            p = subprocess.Popen(exec_args, stdout=subprocess.PIPE)
            st_out = p.communicate()[0].decode("utf-8")
            ret_code = p.returncode

            if ret_code == 0:
                st_parts = st_out.split("MAPQ=")
                mapq_one_res = float(st_parts[1])
                mapq_res[one_range[0]] = mapq_one_res
            else:
                logger.warning('MAPQ returned no data for patient file %s position %s:%s-%s' % (patient_file, chromo,one_range[1], one_range[2]))

        return mapq_res

    @staticmethod
    def execute_freebayes_mp(sargs):
        run_id = sargs[0]
        ref_dir = sargs[1]
        input_dir = sargs[2]
        tmp_dir = sargs[3]
        fasta_ref = sargs[4]
        patient_file = sargs[5]
        chromo = sargs[6]
        pos_ranges = sargs[7]

        count_runs = 0
        for one_range in pos_ranges:
            count_runs = count_runs + 1
            range_arg = '%s:%s-%s'%(chromo, one_range[1], one_range[2])

            out_filename = '%s.%s.%s.vcf'%(run_id,chromo, one_range[0])

            exec_args = ['freebayes', '--region', range_arg, '--gvcf', '--fasta-reference',
                '%s/%s'%(ref_dir,getFastaFileName(fasta_ref)), '%s/%s'%(input_dir, patient_file)]

            logger.info('%d,%s,%s: executing: %s'%(count_runs,chromo, one_range[0],exec_args))

            p = subprocess.Popen(exec_args, stdout=subprocess.PIPE)
            vcf_out = p.communicate()[0].decode("utf-8")
            ret_code = p.returncode
            if ret_code == 0:
                out_file = open('%s%s'%(tmp_dir, out_filename), 'w')
                out_file.write(vcf_out)
                out_file.close()
            else:
                raise RylanProcessException('Failed to execute freebayes on chromo pair. Aborting...')
        pass

    @staticmethod
    def get_position_list_no_flatten(chromo_classif, ref_genome):
        positions_found = []
        for allele in chromo_classif.keys():
            pos_value = chromo_classif[allele][ref_genome]["position"]
            pos_found = (allele,pos_value -1, pos_value)
            positions_found.append(pos_found)

        return positions_found

    @staticmethod
    def get_chromo_positions_info_no_flatten(chromo_classifiers, ref_genome):
        ranges_by_chromo = dict()
        for chromo in chromo_classifiers.keys():
            positions_list = RylanVcfData.get_position_list_no_flatten(chromo_classifiers[chromo], ref_genome)
            ranges_by_chromo[chromo] = positions_list
        return ranges_by_chromo

    @staticmethod
    def get_chromo_positions_info(classifier_rules, ref_genome):
        ranges_by_chromo = {}
        for chromo in classifier_rules.keys():
            positions_list = RylanVcfData.get_position_list(classifier_rules[chromo], ref_genome)
            ranges_by_chromo[chromo] = positions_list
        return ranges_by_chromo

    @staticmethod
    def get_position_list(classifier_rules, ref_genome):
        positions_found = []
        for allele in classifier_rules.keys():
            pos_value = classifier_rules[allele][ref_genome]['position']
            pos_found = (allele.replace('/','_'),pos_value - 1, pos_value)
            positions_found.append(pos_found)

        return positions_found
