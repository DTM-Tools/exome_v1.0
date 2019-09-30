import os
import json
import sys
import traceback
import datetime

from param_stage import *
from vcf_rylan_stage import *

class OutputWriter:

    @staticmethod
    def process_output_stage(pipeline_state):
        rylan_params = pipeline_state['params']
        rylan_vcf_findings = pipeline_state['vcfInfo']
        rylan_annovar_findings = pipeline_state.get('additionalFindings', None)
        rylan_out = OutputWriter(rylan_vcf_findings, rylan_annovar_findings, rylan_params.get_output_directory(),
                                 rylan_params.get_results_file(), rylan_params.get_reference_build(),
                                 rylan_params.get_patient_file(), rylan_params.get_exec_annovar())
        rylan_out.write_output()


    def __init__(self, rylan_container, annovar_container, output_dir, results_file, genref, input_file, exec_annovar):
        self.input_file = input_file
        self.rylan_container = rylan_container
        self.annovar_container = annovar_container
        self.output_dir = output_dir
        self.results_file = results_file
        self.genref = genref
        self.exec_annovar = exec_annovar


    def prepare_findings(self):
        findings_entry = dict()
        delta = datetime.datetime.now() - datetime.datetime(1970, 1, 1)
        ellapsed_seconds = delta.total_seconds()

        findings_entry['runId'] = self.rylan_container.get_run_id()
        findings_entry['timestamp'] = ellapsed_seconds
        findings_entry['fastaRef'] = self.genref
        findings_entry['inputFile'] = self.input_file
        findings_entry['findings'] = self.rylan_container.get_filtered_results()


        if self.annovar_container is not None:
            findings_entry['additionalFindings'] = self.annovar_container.get_annovar_filtered_findings()


        #OutputWriter.add_classification_summary(findings_entry)
        #summary_str = OutputWriter.stringify_summary(findings_entry['summary'])
        #findings_entry['summary_code'] = summary_str

        OutputWriter.pretty_print_findings(findings_entry)

        findings_entry_json = json.dumps(findings_entry)
        return findings_entry_json

    def write_output(self):
        try:

            findings_str = self.prepare_findings()
            if self.results_file != 'stdout':

                out_fname = '%s/%s' % (self.output_dir, self.results_file)
                out_file = open(out_fname, 'w')
                out_file.write(findings_str)
                out_file.close()
            else:
                print(findings_str)

            return
        except Exception as e:
            logger.error('Failed to process rule stage: %s' % str(e))
            traceback.print_exc(file=sys.stdout)
            sys.exit(1)


    '''@staticmethod
    def add_classification_summary(findings_entry):

        findings_summary = findings_entry['findings']

        findings_alleles = sorted(findings_summary.keys())

        summary_dict = dict()

        for found_allele in findings_alleles:
            allele_info = findings_summary[found_allele]

            if allele_info['determination'] == 'classified':

                zygo_code = ''
                if allele_info['zygosity'] == 0:
                    zygo_code = "HOM_REF"
                elif allele_info['zygosity'] == 1:
                    zygo_code = "HET"
                elif allele_info['zygosity'] == 2:
                    zygo_code = 'HOM_VAR'
                else:
                    zygo_code = "ERR"

                allele_entry = dict()

                allele_entry['classification'] = allele_info['classification']
                allele_entry['description'] = allele_info['description']
                allele_entry['zygosity'] = zygo_code
                summary_dict[found_allele] = allele_entry

        findings_entry['summary'] = summary_dict
        return'''

    '''@staticmethod
    def stringify_summary(findings_summary):
        findings_alleles = sorted(findings_summary.keys())
        class_string = ''
        for found_allele in findings_alleles:
            allele_info = findings_summary[found_allele]
            class_string = class_string + found_allele + ':' + allele_info['description'] + ':' + allele_info['classification'] + ':' + allele_info['zygosity'] + ','

        class_string = class_string[:-1]
        return class_string'''


    @staticmethod
    def pretty_print_findings(new_findings):
        print('')
        print('Patient ID: %s\n'%new_findings['inputFile'])
        print('FASTA reference build: %s\n'%new_findings['fastaRef'])
        findings_summary = new_findings['findings']

        findings_alleles = findings_summary.keys()

        for found_allele in findings_alleles:
            allele_info = findings_summary[found_allele]
            if (allele_info['determination'] == "classified"):
                print ("Allele: %s, description: %s, classification: %s"%(found_allele, allele_info['description'], allele_info['classification']) )
            elif (allele_info['determination'] == "filtered"):
                print ("Allele: %s, description: %s, classification: %s FILTERED"%(found_allele, allele_info['description'], allele_info['classification']) )
            else:
                print ("Allele: %s, description: %s, issue: %s"%(found_allele, allele_info['description'], allele_info['determination']))
        print('\n\n\n')
        pass
