

class RylanParams:

    @staticmethod
    def process_params_stage(pipeline_state, rules_file, rules_indel_file, rules_multi_file, patient_file,
                             results_file, ref_build, ref_dir, input_dir, output_dir,
                             cleanup_data = False, num_cores = 4, exec_annovar = False):
        rylan_params = RylanParams(rules_file, rules_indel_file, rules_multi_file, patient_file, results_file, ref_build,
                             ref_dir, input_dir, output_dir, cleanup_data, num_cores, exec_annovar)
        pipeline_state['params'] = rylan_params

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

    def get_rules_file(self):
        return self.rules_file

    def get_rules_indel_file(self):
        return self.rules_indel_file

    def get_rules_multi_file(self):
        return self.rules_multi_file

    def get_patient_file(self):
        return self.patient_file

    def get_results_file(self):
        return self.results_file

    def get_reference_build(self):
        return self.ref_build

    def get_reference_directory(self):
        return self.ref_dir

    def get_input_directory(self):
        return self.input_dir

    def get_output_directory(self):
        return self.output_dir

    def get_cleanup_data(self):
        return self.cleanup_data

    def get_num_cores(self):
        return self.num_cores

    def get_exec_annovar(self):
        return self.exec_annovar
