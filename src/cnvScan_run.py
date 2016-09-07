#!/usr/bin/env python

import sys, os
import pybedtools
import pysam
import vcf
import re
import vcf
import numpy as np
from scipy.stats import mannwhitneyu
import filt_cnvs
import annotate
import argparse
from collections import OrderedDict

class cnv_scan(object):

    def __init__(self, input, output, resources, database):
        self.input = input
        self.output = output
        self.resources = resources
        self.db = database

        self.annotate()
        self.dump()

    def annotate(self):

        self.cnv_anno = {}
        self.cnv_anno, self.cnvs_ordered = filt_cnvs.read_cnvRes(self.input, self.cnv_anno)

        a_cnv = annotate.create_bedTools(self.input)

        # define preprocessors
        p_path = lambda x: os.path.join(self.resources, x)
        p_bed = lambda x: pybedtools.BedTool(p_path(x))
        p_tabix = lambda x: pysam.TabixFile(p_path(x))
        p_vcf = lambda x: vcf.Reader(open(p_path(x)))

        # define annotators
        a_filt_cnvs = lambda x, d: filt_cnvs.db_search(x, d)
        a_gencode = lambda x, d: annotate.gencode_annotate(a_cnv, x, d)
        a_sanger = lambda x, d: annotate.sanger_annotate(a_cnv, x, d)
        a_dgv = lambda x, d: annotate.dgv_annotate(a_cnv, x, d)
        a_dgvFilt2 = lambda x, d: annotate.dgvFilt_annotate(x, d, "DGV_Stringency2")
        a_dgvFilt12 = lambda x, d: annotate.dgvFilt_annotate(x, d, "DGV_Stringency12")
        phastCon = lambda x, d: annotate.phastCon_annotate(x, d)
        a_haplotIdx = lambda x, d: annotate.haploIdx_annotate(x, d)
        a_del1000g = lambda x, d: annotate.del1000g_annotate(x, d)
        a_dup1000g = lambda x, d: annotate.dup1000g_annotate(x, d)
        a_clinVar = lambda x, d: annotate.clinVar_annotate(x, d)
        a_omim = lambda x, d: annotate.omim_annotate(x, d)
        devDisorder = lambda x, d: annotate.devDisorder_annotate(x, d)
        genIntol = lambda x, d: annotate.geneticIntolarance_annotate(x, d)

        # filename, preproc (p), annotate (a)
        queue = [
            (self.db, p_tabix, a_filt_cnvs),
            ("havana_or_ensembl_gencode.v19.annotation.gtf", p_bed, a_gencode),
            ("conrad.et.al.2010_Validated_CNVEs_v5_4Release.tab", p_bed, a_sanger),
            ("dgv_GRCh37_hg19_variants_2014-10-16.tab", p_bed, a_dgv),
            ("cnvMap_stringencyLevel2.bed.gz", p_tabix, a_dgvFilt2),
            ("cnvMap_stringencyLevel12.bed.gz", p_tabix, a_dgvFilt12),
            ("phastConsElements100wayFormatted.bed.gz", p_tabix, phastCon),
            ("haploinsufficiencyindex_withimputation.bed.gz", p_tabix, a_haplotIdx),
            ("union.2010_06.deletions.sites.vcf.gz", p_tabix, a_del1000g),
            ("union.2010_09.TandemDuplications.genotypes.vcf.gz", p_tabix, a_dup1000g),
            ("clinvar_20150106.vcf.gz", p_vcf, a_clinVar),
            ("morbidmap_formatted_onlyHGNC.txt", p_path, a_omim),
            ("cnvScan_DDG2P_freeze_with_gencode19_genomic_coordinates_20141118.txt", p_path, devDisorder),
            ("GeneticIntollarenceScore_RVIS_OERatioPercentile.txt", p_path, genIntol)
        ]

        for (filename, preproc, anno) in queue:
            tmp = preproc(filename)
            self.cnv_anno = anno(tmp, self.cnv_anno)


    def dump(self):

        def gene_names(k, v):
            # 'gene_name': {'SDHDP6':'F', 'RHD':'P', 'C1orf63':'P'} --> "SDHDP6:F|RHD:P|C1orf63:P"
            tmp = v.get('gene_name')
            return '|'.join([k + ':' + v for k, v in tmp.iteritems() ]) if tmp else 'NA'

        def gene_type(k, v):
            # 'gene_type': {'protein_coding': 1, 'processed_transcript': 1} --> "protein_coding;processed_transcript"
            tmp = v.get('gene_type')
            return ";".join(tmp.keys()) if tmp else 'NA'

        def gene_id(k, v):
            tmp = v.get('gene_id')
            return ";".join(tmp.keys()) if tmp else 'NA'

        def exon_count(k, v):
            # 'exon_count': {'ENST00000603639.1': 3, 'ENST00000604864.1': 3} --> "ENST00000603639.1:3|ENST00000604864.1:3"
            tmp = v.get('exon_count')
            return '|'.join([k + ':' + str(v) for k, v in tmp.iteritems()]) if tmp else 'NA'

        annotations = [  # (header, writer)
            ('chr',                     lambda k, v: k.split(":")[0]),
            ('start',                   lambda k, v: k.split(":")[1]),
            ('end',                     lambda k, v: k.split(":")[2]),
            ('cnv_state',               lambda k, v: v.get('CNV_st', 'NA')),
            ('default_score',           lambda k, v: v.get('score', 'NA')),
            ('len',                     lambda k, v : int(k.split(":")[2]) - int(k.split(":")[1])),
            ('inDB_count',              lambda k, v: v.get('inDB_count', 'NA')),
            ('inDBScore_MinMaxMedian',  lambda k, v: v.get('inDB_minmaxmedian', 'NA')),
            ('gene_name',               gene_names),
            ('gene_type',               gene_type),
            ('gene_id',                 gene_id),
            ('exon_count',              exon_count),
            ('UTR',                     lambda k, v: v.get('UTR', 'NA')),
            ('transcript',              lambda k, v: v.get('transcript', 'NA')),
            ('phastConElement_count',   lambda k, v: v.get('phastCon_count', 'NA')),
            ('phastConElement_minMax',  lambda k, v: v.get('phastCon_min_max', 'NA')),
            ('haplo_insufIdx_count',    lambda k, v: v.get('haploIdx_count', 'NA')),
            ('haplo_insufIdx_score',    lambda k, v: v.get('haploIdx_score', 'NA')),
            ('Gene_intolarance_score',  lambda k, v: v.get('GenInTolScore', 'NA')),
            ('sanger_cnv',              lambda k, v: v.get('Sanger_HiRes_CNV', 'NA')),
            ('dgv_cnv',                 lambda k, v: v.get('DGV_CNV', 'NA')),
            ('dgv_varType',             lambda k, v: v.get('DGV_VarType', 'NA')),
            ('dgv_varSubType',          lambda k, v: v.get('DGV_VarSubType', 'NA')),
            ('dgv_pubmedId',            lambda k, v: v.get('DGV_PUBMEDID', 'NA')),
            ('DGV_Stringency2_count',   lambda k, v: v.get('DGV_Stringency2_count', 'NA')),
            ('DGV_Stringency2_PopFreq', lambda k, v: v.get('DGV_Stringency2_popFreq', 'NA')),
            ('DGV_Stringency12_count',  lambda k, v: v.get('DGV_Stringency12_count', 'NA')),
            ('DGV_Stringency12_popFreq',lambda k, v: v.get('DGV_Stringency12_popFreq', 'NA')),
            ('1000g_del',               lambda k, v: v.get('1000G_Del_count', 'NA')),
            ('1000g_ins',               lambda k, v: v.get('1000G_Dup_count', 'NA')),
            ('omim_morbidMap',          lambda k, v: v.get('OMIM', 'NA')),
            ('ddd_mutConsequence',      lambda k, v: v.get('devDis_mutConseq', 'NA')),
            ('ddd_diseaseName',         lambda k, v: v.get('devDis_disName', 'NA')),
            ('ddd_pubmedId',            lambda k, v: v.get('devDis_pubmedID', 'NA')),
            ('clinVar_disease',         lambda k, v: v.get('clindbn', 'NA')),
            ('hgvs_varName',            lambda k, v: v.get('clinhgvs', 'NA')),
            ]

        # re-sorting the dictionary
        cnv_dict = OrderedDict()
        for pos in self.cnvs_ordered:
            cnv_dict[pos] = self.cnv_anno[pos]
            chrom, start, end = pos.split(":")

        with open(self.output, 'w') as out_file:
            header_line = [ header for (header, writer) in annotations ]
            out_file.write("\t".join(header_line)+"\n")

            for key, value in cnv_dict.iteritems():
                line = [ str(writer(key, value)) for (header, writer) in annotations ]
                out_file.write("\t".join( line ) + "\n")


# ============================================================
def main():
        '''
        cnvScan
        '''
        parser = argparse.ArgumentParser(
            description='Annotate CNV prediction with resource data')
        parser.add_argument('-i', '--input', type=str, required=True,
                            help='Input bed file')
        parser.add_argument('-o', '--output', type=str, required=True,
                            help='Output bed file')
        parser.add_argument('-s', '--resources', type=str, required=True,
                            help='Path to resource folder')
        parser.add_argument('-db', '--database', type=str, required=True,
                            help='In-house database file')
        args = parser.parse_args()

        neighbor = cnv_scan(args.input, args.output, args.resources, args.database)


# ============================================================
if __name__ == "__main__":
    sys.exit(main())
