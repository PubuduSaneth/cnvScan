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

        out_file = open(self.output, 'w')

        out_file.write("\t".join(header_line)+"\n")

        for k in self.cnvs_ordered:
            line = k.split(":")
            line.extend([str(self.cnv_anno[k]['CNV_st']), self.cnv_anno[k]['score'], str(int(line[2])-int(line[1])) ])
            line.extend([ str(self.cnv_anno[k]['inDB_count']), str(self.cnv_anno[k]['inDB_minmaxmedian']) ])
            lst_name = []
            exon_c = []

            if self.cnv_anno[k].get('gene_name'):
                for k1 in self.cnv_anno[k]['gene_name']: lst_name.append( ":".join( [k1, self.cnv_anno[k]['gene_name'][k1]] ))
                line.append("|".join(lst_name))
                line.append(";".join(self.cnv_anno[k]['gene_type'].keys()))
                line.append(";".join(self.cnv_anno[k]['gene_id'].keys()))
                if self.cnv_anno[k].get("exon_count"):
                    for k1 in self.cnv_anno[k]['exon_count']: exon_c.append( ":".join( [k1, str(self.cnv_anno[k]['exon_count'][k1])] ))
                    line.append("|".join(exon_c))
                else:
                    line.append("NA")
                if self.cnv_anno[k].get('UTR'):
                    line.append(self.cnv_anno[k]['UTR'])
                else:
                    line.append("NA")
            else:
                self.cnv_anno[k]['gene_name'] = "NA"
                line.append(self.cnv_anno[k]['gene_name'])
                self.cnv_anno[k]['gene_type'] = "NA"
                line.append(self.cnv_anno[k]['gene_type'])
                self.cnv_anno[k]['gene_id'] = "NA"
                line.append(self.cnv_anno[k]['gene_id'])
                self.cnv_anno[k]['exon_count'] = "NA"
                line.append(self.cnv_anno[k]['exon_count'])
                self.cnv_anno[k]['UTR'] = "NA"
                line.append(self.cnv_anno[k]['UTR'])
            if self.cnv_anno[k].get('transcript'):
                line.append(self.cnv_anno[k]['transcript'])
            else:
                self.cnv_anno[k]['transcript'] = "NA"
                line.append(self.cnv_anno[k]['transcript'])
            line.append(str(self.cnv_anno[k]['phastCon_count']))
            line.append(str(self.cnv_anno[k]['phastCon_min_max']))
            line.append(str(self.cnv_anno[k]['haploIdx_count']))
            line.append(str(self.cnv_anno[k]['haploIdx_score']))
            line.append(str(self.cnv_anno[k]['GenInTolScore'])) #
            if self.cnv_anno[k].get('Sanger_HiRes_CNV'):
                line.append(str(self.cnv_anno[k]['Sanger_HiRes_CNV']))
            else:
                self.cnv_anno[k]['Sanger_HiRes_CNV'] = "NA"
                line.append(self.cnv_anno[k]['Sanger_HiRes_CNV'])
            if self.cnv_anno[k].get('DGV_CNV'):
                line.append(str(self.cnv_anno[k]['DGV_CNV']))
                line.append(str(self.cnv_anno[k]['DGV_VarType']))
                line.append(str(self.cnv_anno[k]['DGV_VarSubType']))
                line.append(str(self.cnv_anno[k]['DGV_PUBMEDID']))
            else:
                self.cnv_anno[k]['DGV_CNV'] = "NA"
                self.cnv_anno[k]['DGV_VarType'] = "NA"
                self.cnv_anno[k]['DGV_VarSubType'] = "NA"
                self.cnv_anno[k]['DGV_PUBMEDID'] = "NA"
                line.append(self.cnv_anno[k]['DGV_CNV'])
                line.append(self.cnv_anno[k]['DGV_VarType'])
                line.append(self.cnv_anno[k]['DGV_VarSubType'])
                line.append(self.cnv_anno[k]['DGV_PUBMEDID'])
            if self.cnv_anno[k].get('DGV_Stringency2_count'):
                line.append(str(self.cnv_anno[k]['DGV_Stringency2_count']))   #
                line.append(str(self.cnv_anno[k]['DGV_Stringency2_popFreq'])) #
            else:
                self.cnv_anno[k]['DGV_Stringency2_count'] = "NA"
                self.cnv_anno[k]['DGV_Stringency2_popFreq'] = "NA"
                line.append(str(self.cnv_anno[k]['DGV_Stringency2_count']))
                line.append(str(self.cnv_anno[k]['DGV_Stringency2_popFreq']))
            if self.cnv_anno[k].get('DGV_Stringency12_count'):
                line.append(str(self.cnv_anno[k]['DGV_Stringency12_count']))   #
                line.append(str(self.cnv_anno[k]['DGV_Stringency12_popFreq'])) #
            else:
                self.cnv_anno[k]['DGV_Stringency12_count'] = "NA"
                self.cnv_anno[k]['DGV_Stringency12_popFreq'] = "NA"
                line.append(str(self.cnv_anno[k]['DGV_Stringency12_count']))
                line.append(str(self.cnv_anno[k]['DGV_Stringency12_popFreq']))
            line.append(str(self.cnv_anno[k]['1000G_Del_count']))
            line.append(str(self.cnv_anno[k]['1000G_Dup_count']))
            line.append(self.cnv_anno[k]['OMIM'])
            line.append(self.cnv_anno[k]['devDis_mutConseq'])
            line.append(self.cnv_anno[k]['devDis_disName'])
            line.append(self.cnv_anno[k]['devDis_pubmedID'])
            line.append(str(self.cnv_anno[k]['clindbn']))
            line.append(str(self.cnv_anno[k]['clinhgvs']))
            out_file.write("\t".join(line) + "\n")
            #print "\t".join(line), len(line)


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
