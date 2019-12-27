import sys, getopt
import datetime
from pipeline import *
import os
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger(__name__)

def print_usage():
    print('usage: main.py [-a] [-n] -m <num_cores> -r <refdir> -o <outputdir> -i <inputdir> -p <patientfile.dedup.bam> -g <hg19|grch38>')


def main(argv):
    outputdir = ''
    refdir = ''
    patientfile = ''
    inputdir = ''
    genref = 'hg19'
    cleanupdata = True
    numcores = 1
    run_annovar = False
    try:
        opts, args = getopt.getopt(argv,"anho:i:r:p:m:g:d:",["outputdir=", "inputdir=", "refdir=", "patientfile=", "mcores=", "genref="])
    except getopt.GetoptError:
        print_usage()
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt == '-n':
            cleanupdata = False
        elif opt == '-a':
            run_annovar = True
        elif opt in ("-r", "--refdir"):
            refdir = arg
        elif opt in ("-i", "--inputdir"):
            inputdir = arg
        elif opt in ("-p", "--patientfile"):
            patientfile = arg
        elif opt in ("-o", "--outputdir"):
            outputdir = arg
        elif opt in ("-g", "--genref"):
            genref = arg
        elif opt in ("-m", "--mcores"):
            numcores = int(arg)
        else:
            #print(opt)
            print('Unrecognized option.')
            print_usage()
            sys.exit(0)

    if run_annovar is True:
        run_annovar_str = 'True'
    else:
        run_annovar_str = 'False'

    logger.info('FASTA ref dir is %s'%(refdir))
    logger.info('FASTA reference is %s'%(genref))
    logger.info('Output dir is: %s'%(outputdir))
    logger.info('Patient file is %s'%(patientfile))
    logger.info('Input dir is %s'%(inputdir))
    logger.info('Clean-up temporary VCFs: %s'%(cleanupdata))
    logger.info('Number of CPU cores: %d'%(numcores))
    logger.info('Run Annovar: %s'%(run_annovar_str))

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
