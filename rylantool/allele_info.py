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
c_hg19_pos = c_chromo + 1
c_gr38_pos = c_hg19_pos + 1
c_A_cat = c_gr38_pos + 1
c_C_cat = c_A_cat + 1
c_G_cat = c_C_cat + 1
c_T_cat = c_G_cat + 1
c_A_class= c_T_cat + 1
c_C_class = c_A_class + 1
c_G_class =  c_C_class + 1
c_T_class = c_G_class + 1
c_DP_REF = c_T_class + 1
c_MAPQ_REF = c_DP_REF + 1
c_QUAL = c_MAPQ_REF + 1
c_DP = c_QUAL + 1
c_AO = c_DP + 1
c_AO_DP = c_AO + 1
c_Cell = c_AO_DP + 1
c_Type = c_Cell + 1
c_Gensubtype =c_Type + 1
c_Asubtype = c_Gensubtype + 1
c_Csubtype = c_Asubtype + 1
c_Gsubtype = c_Csubtype + 1
c_Tsubtype = c_Gsubtype + 1

def assemble_filter(line_data):
    filter_info = {}
    filter_info['allele'] = line_data[c_allele]
    filter_info['description'] = line_data[c_desc]

    fr_col = []
    fdp_info = {}
    fdp_info['filterId'] = 'DP'
    fdp_info['criteria'] = line_data[c_DP_REF].lower()
    fdp_info['action'] = 'flag'
    fr_col.append(fdp_info)

    fmq_info = {}
    fmq_info['filterId'] = 'MAPQ'
    fmq_info['criteria'] = line_data[c_MAPQ_REF].lower()
    fmq_info['action'] = 'flag'
    fr_col.append(fmq_info)

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

    fv_col.append(fmq_info)

    filter_info['variant'] = fv_col
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

    if line_data[c_A_cat] == "ref":
        reference = ('A', line_data[c_A_class])
    elif line_data[c_A_cat] == "var":
        var_data = ('A', line_data[c_A_class])
        variants.append(var_data)
    elif line_data[c_A_cat] == "abn":
        abnormalities.append('A')

    if line_data[c_C_cat] == "ref":
        reference = ('C', line_data[c_C_class])
    elif line_data[c_C_cat] == "var":
        var_data = ('C', line_data[c_C_class])
        variants.append(var_data)
    elif line_data[c_C_cat] == "abn":
        abnormalities.append('C')

    if line_data[c_G_cat] == "ref":
        reference = ('G', line_data[c_G_class])
    elif line_data[c_G_cat] == "var":
        var_data = ('G', line_data[c_G_class])
        variants.append(var_data)
    elif line_data[c_G_cat] == "abn":
        abnormalities.append('G')

    if line_data[c_T_cat] == "ref":
        reference = ('T', line_data[c_T_class])
    elif line_data[c_T_cat] == "var":
        var_data = ('T', line_data[c_T_class])
        variants.append(var_data)
    elif line_data[c_T_cat] == "abn":
        abnormalities.append('T')

    hg19_info['position'] = int(line_data[c_hg19_pos])
    grch38_info['position'] = int(line_data[c_gr38_pos])

    ref_struct = {}
    ref_struct['nucleotide'] = reference[0]
    ref_struct['classification'] = reference[1]
    ref_struct['action'] = 'accept'

    hg19_info['reference'] = ref_struct
    grch38_info['reference'] = ref_struct

    var_struct = []

    for var_info in variants:
        var_sub_struct = {}
        var_sub_struct['nucleotide'] = var_info[0]
        var_sub_struct['classification'] = var_info[1]
        var_sub_struct['action'] = 'accept'
        var_struct.append(var_sub_struct)

    hg19_info['variants'] = var_struct
    grch38_info['variants'] = var_struct

    abn_struct = []

    for abn_info in abnormalities:
        abn_sub_struct = {}
        abn_sub_struct['nucleotide'] = abn_info
        abn_sub_struct['action'] = 'flag'
        abn_struct.append(abn_sub_struct)

    hg19_info['abnormalities'] = abn_struct
    grch38_info['abnormalities'] = abn_struct

    hg19_spdi = list()
    grch38_spdi = list()
    for a_variant in variants:
        spdi_code_hg19 = calculate_spdi_code(info['chromosome'], hg19_info['position'], reference[0], a_variant[0])
        hg19_spdi.append(spdi_code_hg19)
        spdi_code_grch38 = calculate_spdi_code(info['chromosome'], grch38_info['position'], reference[0], a_variant[0])
        grch38_spdi.append(spdi_code_grch38)

    hg19_info['spdi'] = hg19_spdi
    grch38_info['spdi'] = grch38_spdi

    info['hg19'] = hg19_info
    info['grch38'] = grch38_info

    info['cell'] = line_data[c_Cell]
    info['type'] = line_data[c_Type]
    info['genSubtype'] = line_data[c_Gensubtype]
    info['aSubtype'] = line_data[c_Asubtype]
    info['cSubtype'] = line_data[c_Csubtype]
    info['gSubtype'] = line_data[c_Gsubtype]
    info['tSubtype'] = line_data[c_Tsubtype]

    return info


def read_rylan_allele_info(filein):

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
            line_data = lines[i].replace('\n','').split(",")
            logger.info('Processed data: %s' % line_data)
            entry_data = assemble_info(line_data)
            filter_data = assemble_filter(line_data)
            entry_collection.append(entry_data)
            filter_collection.append(filter_data)


        logger.info('Done processing allele data.')
        return entry_collection, filter_collection

    except Exception as e:
        logger.error(e)
        raise RylanProcessException('Failed to process allele data.')
