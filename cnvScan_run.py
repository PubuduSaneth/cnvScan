#!/usr/bin/env python

#!/usr/bin/env python

import sys
import pybedtools
import pysam
import vcf
import re
import vcf
import numpy as np
from scipy.stats import mannwhitneyu
import filt_cnvs
import annotate


resource_dir = "/home/saneth/Documents/cnvFilt_proj/resources"
cnv_res_file = sys.argv[1]


db_file = pysam.TabixFile(sys.argv[3])

cnv_anno = {}
cnv_anno, cnvs_ordered = filt_cnvs.read_cnvRes(cnv_res_file, cnv_anno)
cnv_anno = filt_cnvs.db_search(db_file, cnv_anno) # cnv_anno = filt_cnvs.db_search(db_file, db_id, cnv_anno) #


a_cnv = annotate.create_bedTools(sys.argv[1])
b_gencode = pybedtools.BedTool(resource_dir+"/gencode19/havana_or_ensembl_gencode.v19.annotation.gtf")
c_conradCNV = pybedtools.BedTool(resource_dir+"/conrad.et.al.2010_Validated_CNVEs_v5_4Release.tab")
d_dgvCNV = pybedtools.BedTool(resource_dir+"/dgv_GRCh37_hg19_variants_2014-10-16.tab")
d_dgvFiltsCNV_l2 = pysam.TabixFile(resource_dir+"/dgv-filtered/cnvMap_stringencyLevel2.bed.gz")
d_dgvFiltsCNV_l12 = pysam.TabixFile(resource_dir+"/dgv-filtered/cnvMap_stringencyLevel12.bed.gz")
e_phastCon = pysam.TabixFile(resource_dir+"/PhastCon/phastConsElements100wayFormatted.bed.gz")
f_haploIdx = pysam.TabixFile(resource_dir+"/haploinsufficiencyindex/haploinsufficiencyindex_withimputation.bed.gz")
g_del1000g_delFile = pysam.TabixFile(resource_dir+"/1000GSVs/1000GCNV/union.2010_06.deletions.sites.vcf.gz")
h_dup1000g_delFile = pysam.TabixFile(resource_dir+"/1000GSVs/1000GCNV/union.2010_09.TandemDuplications.genotypes.vcf.gz")
i_clinVar_reader = vcf.Reader(open('/home/saneth/Documents/cnvFilt_proj/resources/clinvar_20150106.vcf.gz', 'r'))
j_omim_file = resource_dir+"/OMIM/morbidmap_formatted_onlyHGNC.txt"
h_devDis_file = resource_dir+"/ddg2p_20141118/cnvScan_DDG2P_freeze_with_gencode19_genomic_coordinates_20141118.txt"
i_genIntol_file = resource_dir+"/GeneticIntollarenceScore/GeneticIntollarenceScore_RVIS_OERatioPercentile.txt"

cnv_anno = annotate.gencode_annotate(a_cnv, b_gencode, cnv_anno)
cnv_anno = annotate.sanger_annotate(a_cnv, c_conradCNV, cnv_anno)
cnv_anno = annotate.dgv_annotate(a_cnv, d_dgvCNV, cnv_anno)
cnv_anno = annotate.dgvFilt_annotate(d_dgvFiltsCNV_l2, cnv_anno, "DGV_Stringency2")
cnv_anno = annotate.dgvFilt_annotate(d_dgvFiltsCNV_l12, cnv_anno, "DGV_Stringency12")
cnv_anno = annotate.phastCon_annotate(e_phastCon, cnv_anno)
cnv_anno = annotate.haploIdx_annotate(f_haploIdx, cnv_anno)
cnv_anno = annotate.geneticIntolarance_annotate(i_genIntol_file, cnv_anno)
cnv_anno = annotate.del1000g_annotate(g_del1000g_delFile, cnv_anno)
cnv_anno = annotate.dup1000g_annotate(h_dup1000g_delFile, cnv_anno)
cnv_anno = annotate.clinVar_annotate(i_clinVar_reader, cnv_anno)
cnv_anno = annotate.omim_annotate(j_omim_file, cnv_anno)
cnv_anno = annotate.devDisorder_annotate(h_devDis_file, cnv_anno)

header_line= ["chr", "start", "end", "cnv_state", "score","len"]
header_line.extend(["inDB_count", "inDB_MinMaxMedian"])
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

out_file = open(sys.argv[2], 'w')

out_file.write("\t".join(header_line)+"\n")

for k in cnvs_ordered:
    line = k.split(":")
    line.extend([str(cnv_anno[k]['CNV_st']), cnv_anno[k]['score'], str(int(line[2])-int(line[1])) ])
    line.extend([ str(cnv_anno[k]['inDB_count']), str(cnv_anno[k]['inDB_minmaxmedian']) ])
    lst_name = []
    exon_c = []

    if cnv_anno[k].get('gene_name'):
        for k1 in cnv_anno[k]['gene_name']: lst_name.append( ":".join( [k1, cnv_anno[k]['gene_name'][k1]] ))
        line.append("|".join(lst_name))
        line.append(";".join(cnv_anno[k]['gene_type'].keys()))
        line.append(";".join(cnv_anno[k]['gene_id'].keys()))
        if cnv_anno[k].get("exon_count"):
            for k1 in cnv_anno[k]['exon_count']: exon_c.append( ":".join( [k1, str(cnv_anno[k]['exon_count'][k1])] ))
            line.append("|".join(exon_c))
        else:
            line.append("NA")
        if cnv_anno[k].get('UTR'):
            line.append(cnv_anno[k]['UTR'])
        else:
            line.append("NA")
    else:
        cnv_anno[k]['gene_name'] = "NA"
        line.append(cnv_anno[k]['gene_name'])
        cnv_anno[k]['gene_type'] = "NA"
        line.append(cnv_anno[k]['gene_type'])
        cnv_anno[k]['gene_id'] = "NA"
        line.append(cnv_anno[k]['gene_id'])
        cnv_anno[k]['exon_count'] = "NA"
        line.append(cnv_anno[k]['exon_count'])
        cnv_anno[k]['UTR'] = "NA"
        line.append(cnv_anno[k]['UTR'])
    if cnv_anno[k].get('transcript'):
        line.append(cnv_anno[k]['transcript'])
    else:
        cnv_anno[k]['transcript'] = "NA"
        line.append(cnv_anno[k]['transcript'])
    line.append(str(cnv_anno[k]['phastCon_count']))
    line.append(str(cnv_anno[k]['phastCon_min_max']))
    line.append(str(cnv_anno[k]['haploIdx_count']))
    line.append(str(cnv_anno[k]['haploIdx_score']))
    line.append(str(cnv_anno[k]['GenInTolScore'])) #
    if cnv_anno[k].get('Sanger_HiRes_CNV'):
        line.append(str(cnv_anno[k]['Sanger_HiRes_CNV']))
    else:
        cnv_anno[k]['Sanger_HiRes_CNV'] = "NA"
        line.append(cnv_anno[k]['Sanger_HiRes_CNV'])
    if cnv_anno[k].get('DGV_CNV'):
        line.append(str(cnv_anno[k]['DGV_CNV']))
        line.append(str(cnv_anno[k]['DGV_VarType']))
        line.append(str(cnv_anno[k]['DGV_VarSubType']))
        line.append(str(cnv_anno[k]['DGV_PUBMEDID']))
    else:
        cnv_anno[k]['DGV_CNV'] = "NA"
        cnv_anno[k]['DGV_VarType'] = "NA"
        cnv_anno[k]['DGV_VarSubType'] = "NA"
        cnv_anno[k]['DGV_PUBMEDID'] = "NA"
        line.append(cnv_anno[k]['DGV_CNV'])
        line.append(cnv_anno[k]['DGV_VarType'])
        line.append(cnv_anno[k]['DGV_VarSubType'])
        line.append(cnv_anno[k]['DGV_PUBMEDID'])
    if cnv_anno[k].get('DGV_Stringency2_count'):
        line.append(str(cnv_anno[k]['DGV_Stringency2_count']))   #
        line.append(str(cnv_anno[k]['DGV_Stringency2_popFreq'])) #
    else:
        cnv_anno[k]['DGV_Stringency2_count'] = "NA"
        cnv_anno[k]['DGV_Stringency2_popFreq'] = "NA"
        line.append(str(cnv_anno[k]['DGV_Stringency2_count']))
        line.append(str(cnv_anno[k]['DGV_Stringency2_popFreq']))
    if cnv_anno[k].get('DGV_Stringency12_count'):
        line.append(str(cnv_anno[k]['DGV_Stringency12_count']))   #
        line.append(str(cnv_anno[k]['DGV_Stringency12_popFreq'])) #
    else:
        cnv_anno[k]['DGV_Stringency12_count'] = "NA"
        cnv_anno[k]['DGV_Stringency12_popFreq'] = "NA"
        line.append(str(cnv_anno[k]['DGV_Stringency12_count']))
        line.append(str(cnv_anno[k]['DGV_Stringency12_popFreq']))    
    line.append(str(cnv_anno[k]['1000G_Del_count']))
    line.append(str(cnv_anno[k]['1000G_Dup_count']))
    line.append(cnv_anno[k]['OMIM'])    
    line.append(cnv_anno[k]['devDis_mutConseq'])
    line.append(cnv_anno[k]['devDis_disName'])
    line.append(cnv_anno[k]['devDis_pubmedID'])
    line.append(str(cnv_anno[k]['clindbn']))
    line.append(str(cnv_anno[k]['clinhgvs']))
    out_file.write("\t".join(line) + "\n")
    #print "\t".join(line), len(line)
