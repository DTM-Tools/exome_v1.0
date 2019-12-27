import sys
import os
import logging
import json

from rylanexceptions import *
from utils import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

c_allele = 0
c_desc = c_allele + 1
c_enabled = c_desc + 1
c_system = c_enabled + 1
c_gene = c_system + 1
c_chromo = c_gene + 1
c_hg19_pos_st = c_chromo + 1
c_hg19_pos_end = c_hg19_pos_st + 1
c_gr38_pos_st = c_hg19_pos_end + 1
c_gr38_pos_end = c_gr38_pos_st + 1
c_Ref = c_gr38_pos_end + 1
c_Alt = c_Ref + 1
c_Ref_class= c_Alt + 1
c_Alt_class = c_Ref_class + 1
c_Type =  c_Alt_class + 1
c_DP_Ref = c_Type + 1
c_MAPQ = c_DP_Ref + 1
c_MQM = c_MAPQ + 1
c_MQMR = c_MQM + 1
c_QUAL = c_MQMR + 1
c_DP = c_QUAL + 1
c_AO = c_DP + 1
c_AO_DP = c_AO + 1
c_QR = c_AO_DP + 1
c_QA = c_QR + 1
c_Cell = c_QA + 1
c_Function = c_Cell + 1
c_Gensubtype = c_Function + 1
c_Refsubtype = c_Gensubtype + 1
c_Altsubtype = c_Refsubtype + 1

def assemble_filter(line_data):
    filter_info = {}
    filter_info['allele'] = line_data[c_allele]
    filter_info['description'] = line_data[c_desc]

    fr_col = []
    fdp_info = {}
    fdp_info['filterId'] = 'DP'
    fdp_info['criteria'] = line_data[c_DP_Ref].lower()
    fdp_info['action'] = 'flag'
    fr_col.append(fdp_info)

    fmqm_info = {}
    fmqm_info['filterId'] = 'MQM'
    fmqm_info['criteria'] = line_data[c_MQM].lower()
    fmqm_info['action'] = 'flag'
    #fr_col.append(fmqm_info)

    fmqmr_info = {}
    fmqmr_info['filterId'] = 'MQMR'
    fmqmr_info['criteria'] = line_data[c_MQMR].lower()
    fmqmr_info['action'] = 'flag'
    #fr_col.append(fmqmr_info)

    fqr_info = {}
    fqr_info['filterId'] = 'QR'
    fqr_info['criteria'] = line_data[c_QR].lower()
    fqr_info['action'] = 'flag'
    fr_col.append(fqr_info)

    fqa_info = {}
    fqa_info['filterId'] = 'QA'
    fqa_info['criteria'] = line_data[c_QA].lower()
    fqa_info['action'] = 'flag'
    #fr_col.append(fqa_info)

    filter_info['reference'] = fr_col

    fv_col = []

    vq_info = {}
    vq_info['filterId'] = 'QUAL'
    vq_info['criteria'] = line_data[c_QUAL].lower()
    vq_info['action'] = 'flag'
    fv_col.append(vq_info)

    vdp_info = {}
    vdp_info['filterId'] = 'DP'
    vdp_info['criteria'] = line_data[c_DP].lower()
    vdp_info['action'] = 'flag'
    fv_col.append(vdp_info)

    vao_info = {}
    vao_info['filterId'] = 'AO'
    vao_info['criteria'] = line_data[c_AO].lower()
    vao_info['action'] = 'flag'
    fv_col.append(vao_info)

    vao_dp_info = {}
    vao_dp_info['filterId'] = 'AO/DP'
    vao_dp_info['criteria'] = line_data[c_AO_DP].lower()
    vao_dp_info['action'] = 'flag'
    fv_col.append(vao_dp_info)

    fv_col.append(fmqm_info)
    fv_col.append(fmqmr_info)
    fv_col.append(fqr_info)
    fv_col.append(fqa_info)

    filter_info['variant'] = fv_col

    fhv_col = []

    fhv_col.append(vq_info)
    fhv_col.append(vdp_info)
    fhv_col.append(vao_info)
    fhv_col.append(vao_dp_info)
    fhv_col.append(fqa_info)

    filter_info['homozygous_variant'] = fhv_col

    return filter_info


def assemble_info(line_data):
    info = {}
    info['allele'] = line_data[c_allele]
    info['description'] = line_data[c_desc]
    info['system'] = line_data[c_system]
    info['gene'] = line_data[c_gene]
    info['chromosome'] = line_data[c_chromo]
    info['enabled'] = bool(line_data[c_enabled])
    hg19_info = {}
    grch38_info = {}

    reference = ()
    variants= []
    abnormalities = []


    reference = line_data[c_Ref]
    alternate = line_data[c_Alt]
    indel_type = line_data[c_Type]

    hg19_info['position_start'] = int(line_data[c_hg19_pos_st])
    hg19_info['position_end'] = int(line_data[c_hg19_pos_end])
    grch38_info['position_start'] = int(line_data[c_gr38_pos_st])
    grch38_info['position_end'] = int(line_data[c_gr38_pos_end])

    ref_struct = {}
    ref_struct['nucleotides'] = reference
    ref_struct['class'] = line_data[c_Ref_class]

    alt_struct = {}
    alt_struct['nucleotides'] = alternate
    alt_struct['class'] = line_data[c_Alt_class]

    hg19_info['reference'] = ref_struct
    grch38_info['reference'] = ref_struct

    hg19_info['alternate'] = alt_struct
    grch38_info['alternate'] = alt_struct

    if indel_type == "ins":
        spdi_code_hg19 = calculate_spdi_code_ins(info['chromosome'], hg19_info['position_end'], reference, alternate)
        spdi_code_grch38 = calculate_spdi_code_ins(info['chromosome'], grch38_info['position_end'], reference, alternate)

    elif indel_type == "del":
        spdi_code_hg19 = calculate_spdi_code_del(info['chromosome'], hg19_info['position_start'], reference, alternate)
        spdi_code_grch38 = calculate_spdi_code_del(info['chromosome'], grch38_info['position_start'], reference, alternate)

    hg19_info['spdi'] = [spdi_code_hg19]
    grch38_info['spdi'] = [spdi_code_grch38]

    info['hg19'] = hg19_info
    info['grch38'] = grch38_info

    info['cell'] = line_data[c_Cell]
    info['function'] = line_data[c_Function]
    info['genSubtype'] = line_data[c_Gensubtype]
    info['refSubtype'] = line_data[c_Refsubtype]
    info['altSubtype'] = line_data[c_Altsubtype]

    return info

def read_rylan_indel_info(filein):

    try:
        logger.info('Reading CSV file...')
        fin = open(filein)

        lines = []
        for line in fin:
            lines.append(line)

        fin.close()

        entry_collection = []
        filter_collection = []

        for i in range(1,len(lines)):
            line_data = lines[i].split(",")
            logger.info('Processed data: %s' % line_data)
            entry_data = assemble_info(line_data)
            filter_data = assemble_filter(line_data)
            entry_collection.append(entry_data)
            filter_collection.append(filter_data)

        logger.info('Done processing indel data.')

        return entry_collection, filter_collection

    except Exception as e:
        logger.error(e)
        raise RylanProcessException('Failed to process indel data.')
