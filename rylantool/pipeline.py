
import logging
import sys
import json
import traceback

from param_stage import *
from rules_stage import *
from vcf_rylan_stage import *
from vcf_rylan_indel_stage import *
from output_stage import *
from rylanexceptions import *
from phasing_stage import *
from additional_findings_stage import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class RylanPipeline:

    def __init__(self, rules_file, rules_indel_file, rules_multi_file, patient_file, results_file, ref_build,
                 ref_dir, input_dir, output_dir, cleanup_data = False, num_cores = 1, exec_annovar = False):
        self.rules_file = rules_file
        self.rules_indel_file = rules_indel_file
        self.rules_multi_file = rules_multi_file
        self.patient_file = patient_file
        self.results_file = results_file
        self.ref_build = ref_build
        self.ref_dir = ref_dir
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.cleanup_data = cleanup_data
        self.num_cores = num_cores
        self.exec_annovar = exec_annovar

    def execute_pipeline(self):
        try:
            logger.info('Starting pipeline execution...')

            self.pipeline_state = dict()

            # Step 0. Process parameters...
            print('Ref dir: %s' % self.ref_dir)
            RylanParams.process_params_stage(self.pipeline_state,
                                            self.rules_file, self.rules_indel_file, self.rules_multi_file, self.patient_file,
                                            self.results_file, self.ref_build, self.ref_dir,
                                            self.input_dir, self.output_dir, self.cleanup_data,
                                            self.num_cores, self.exec_annovar)

            # Step 1. Load the rules...
            RylanRules.process_rules_stage(self.pipeline_state)

            # Step 2. Generate VCF Files...
            RylanVcfData.process_vcf_stage(self.pipeline_state)

            # Step 3. Process indels...
            RylanInDelVcfData.process_vcf_indel_stage(self.pipeline_state)

            #Step 4. Process multiples...
            RylanPhasing.process_phasing_stage(self.pipeline_state)

            if self.exec_annovar == True:
                #Step 5. Process anovar...
                RylanAdditionalFindings.process_additional_findings_stage(self.pipeline_state)

            #Last step: output findings
            OutputWriter.process_output_stage(self.pipeline_state)

            # Clean-up
            if self.cleanup_data == True:
                self.pipeline_state['vcfInfo'].cleanup_vcf_files()

            print(self.pipeline_state)
            logger.info('Finished executing RyLAN pipeline.')

        except Exception as e:
            logger.error(str(e))
            traceback.print_exc(file=sys.stdout)
            raise RylanPipelineException(str(e))
