import sys, getopt
import datetime
from pipeline import *
import os
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger(__name__)


def main(argv):
    inputdir = os.environ.get('INPUT_DIR')
    outputdir = os.environ.get('OUTPUT_DIR')
    refdir = os.environ.get('REF_DIR')
    patientfile = os.environ.get('PATIENT_FILE')
    genref =  os.environ.get('REF_BUILD','hg19')
    cleanupdata = True
    numcores = int(os.environ.get('NUMCORES','2'))
    run_annovar_str = os.environ.get('RUN_ANNOVAR','True')

    logger.info('FASTA ref dir is %s'%(refdir))
    logger.info('FASTA reference is %s'%(genref))
    logger.info('Output dir is: %s'%(outputdir))
    logger.info('Patient file is %s'%(patientfile))
    logger.info('Input dir is %s'%(inputdir))
    logger.info('Clean-up temporary VCFs: %s'%(cleanupdata))
    logger.info('Number of CPU cores: %d'%(numcores))
    logger.info('Run Annovar: %s'%(run_annovar_str))

    if run_annovar_str.lower() == 'true':
        run_annovar = True
    else:
        run_annovar = False

    rylan_rules_file = 'ChromoList.csv'
    rylan_indel_file = 'ChromoInDelList.csv'
    rules_multi_file = 'Multi.csv'
    results_file = '%s.rylanout.json' % patientfile

    rylan_pipeline = RylanPipeline(rylan_rules_file, rylan_indel_file, rules_multi_file, patientfile,
                                   results_file, genref, refdir, inputdir, outputdir,
                                   cleanup_data=cleanupdata, num_cores=numcores, exec_annovar=run_annovar)
    logger.info('Executing pipeline...')
    rylan_pipeline.execute_pipeline()

    logger.info('Done processing RyLAN job.')
    sys.exit(0)

if __name__ == "__main__":
   main(sys.argv[1:])
