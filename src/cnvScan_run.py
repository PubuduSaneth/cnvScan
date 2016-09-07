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
        header_line= ["chr", "start", "end", "cnv_state", "default_score","len"]
        header_line.extend(["inDB_count", "inDBScore_MinMaxMedian"])
        header_line.extend(["gene_name", "gene_type", "gene_id", "exon_count", "UTR", "transcript"])
        header_line.extend(["phastConElement_count", "phastConElement_minMax"])
        header_line.extend(["haplo_insufIdx_count", "haplo_insufIdx_score"])
        header_line.append("Gene_intolarance_score")
        header_line.append("sanger_cnv")
        header_line.extend(["dgv_cnv", "dgv_varType", "dgv_varSubType", "dgv_pubmedId", 'DGV_Stringency2_count', 'DGV_Stringency2_PopFreq', 'DGV_Stringency12_count', 'DGV_Stringency12_popFreq'])
        header_line.extend(["1000g_del","1000g_ins"])
        header_line.append("omim_morbidMap")
        header_line.extend(["ddd_mutConsequence", "ddd_diseaseName", "ddd_pubmedId"])
        header_line.extend(["clinVar_disease", "hgvs_varName"])

        with open(self.output, 'w') as out_file:

            out_file.write("\t".join(header_line)+"\n")

            # re-sorting the dictionary
            cnv_dict = OrderedDict()
            for pos in self.cnvs_ordered:
                cnv_dict[pos] = self.cnv_anno[pos]
                chrom, start, end = pos.split(":")

            for key, value in cnv_dict.iteritems():

                line = []
                chrom, start, end = key.split(":")
                line += [chrom]
                line += [start]
                line += [end]
                line += [value.get('CNV_st', 'NA')]
                line += [value.get('score', 'NA')]
                line += [ int(end) - int(start) ]
                line += [value.get('inDB_count', 'NA')]
                line += [value.get('inDB_minmaxmedian', 'NA')]

                # for all: if key does exist, write "NA"
                # 'gene_name': {'SDHDP6':'F', 'RHD':'P', 'C1orf63':'P'} --> "SDHDP6:F|RHD:P|C1orf63:P"
                gene_names = value.get('gene_name')
                line += ['|'.join([k + ':' + v for k, v in gene_names.iteritems() ]) if gene_names else 'NA']

                # 'gene_type': {'protein_coding': 1, 'processed_transcript': 1} --> "protein_coding;processed_transcript"
                gene_type = value.get('gene_type')
                line += [";".join(gene_type.keys()) if gene_type else 'NA']

                gene_id = value.get('gene_id')
                line += [";".join(gene_id.keys()) if gene_id else 'NA']

                # 'exon_count': {'ENST00000603639.1': 3, 'ENST00000604864.1': 3} --> "ENST00000603639.1:3|ENST00000604864.1:3"
                exon_count = value.get('exon_count')
                line += ['|'.join([k + ':' + str(v) for k, v in exon_count.iteritems()]) if exon_count else 'NA']

                line += [value.get('UTR', 'NA')]
                line += [value.get('transcript', 'NA')]
                line += [value.get('phastCon_count', 'NA')]
                line += [value.get('phastCon_min_max', 'NA')]
                line += [value.get('haploIdx_count', 'NA')]
                line += [value.get('haploIdx_score', 'NA')]
                line += [value.get('GenInTolScore', 'NA')]
                line += [value.get('Sanger_HiRes_CNV', 'NA')]
                line += [value.get('DGV_CNV', 'NA')]
                line += [value.get('DGV_VarType', 'NA')]
                line += [value.get('DGV_VarSubType', 'NA')]
                line += [value.get('DGV_PUBMEDID', 'NA')]
                line += [value.get('DGV_Stringency2_count', 'NA')]
                line += [value.get('DGV_Stringency2_popFreq', 'NA')]
                line += [value.get('DGV_Stringency12_count', 'NA')]
                line += [value.get('DGV_Stringency12_popFreq', 'NA')]
                line += [value.get('1000G_Del_count', 'NA')]
                line += [value.get('1000G_Dup_count', 'NA')]
                line += [value.get('OMIM', 'NA')]
                line += [value.get('devDis_mutConseq', 'NA')]
                line += [value.get('devDis_disName', 'NA')]
                line += [value.get('devDis_pubmedID', 'NA')]
                line += [value.get('clindbn', 'NA')]
                line += [value.get('clinhgvs', 'NA')]

                out_file.write("\t".join( [str(elem) for elem in line] ) + "\n")


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
