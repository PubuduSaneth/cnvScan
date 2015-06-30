#!/usr/bin/env python

import re
import pybedtools
import vcf
import numpy as np
from scipy.stats import mannwhitneyu
import pysam

def read_cnvRes(cnv_res_file, cnv_anno):
    cnvs_ordered = []
    
    with open(cnv_res_file, 'r') as f:
        #next(f)
        for line in f:
            l = line.rstrip("\n").split("\t")
            l[0] = l[0] #l[0] = l[0][3:]
            key = ":".join(l[0:3])
            cnv_anno[key] = {}
            cnv_anno[key]['len'] = int(l[2]) - int(l[1])
            if l[3] == "deletion":
                l[3] = 1
            elif l[3] == "duplication":
                l[3] = 3
            cnv_anno[key]['CNV_st'] = int(l[3])
            cnv_anno[key]['score'] = l[4]
            cnvs_ordered.append(key)
    return cnv_anno, cnvs_ordered
    pass

def het_hom_count_2(cnv_anno, vcf_file):
    for key in cnv_anno.keys():
        if cnv_anno[key]['CNV_st'] < 2:  #== "deletion":
            het_count = 0
            hom_count = 0
            line = key.split(":")
            r=[]
            #print "%r" %line
            try:
                for row in vcf_file.fetch(str(line[0]), int(line[1]), int(line[2])):
                    r = row.split("\t") 
                    if r[6] == "PASS":
                        if re.search(r'het',r[7]):
                            het_count = het_count +1
                            #print "HET::: %r" %r
                        elif re.search(r'hom',r[7]):
                            hom_count = hom_count + 1
                            #print "HOM::: %r" %r
                if hom_count > 0:
                    hh_ratio = round(float(het_count) / hom_count, 2)
                else:
                    hh_ratio = "NA"
            except ValueError:
                het_count = "NA"
                hom_count = "NA"
                hh_ratio = "NA"
        else:
            het_count = "NA"
            hom_count = "NA"
            hh_ratio = "NA"
            
        cnv_anno[key]['Het_count'] = het_count
        cnv_anno[key]['Hom_count'] = hom_count
        cnv_anno[key]['Het/Hom_ratio'] = hh_ratio
    return cnv_anno    
    pass


def het_home_count(cnv_anno, vcf_file):
    
    vcf_reader1 = vcf.Reader(open(vcf_file, 'r'))

    for key in cnv_anno.keys():
        line = key.split(":")
        het_count = 0
        hom_count = 0
        if cnv_anno[key]['CNV_st'] == "deletion":
            #print "***** %r %r %r" %(line[0], line[1], line[2])
            try:
                for record in vcf_reader1.fetch(line[0], int(line[1]), int(line[2])):
                    #print "@@@@", record
                    if record.num_het != 0 and len(record.FILTER) == 0:
                    #if record.num_het != 0:
                        het_count = het_count +1
                        #print "###", het_count
                    else:
                        hom_count += 1
                if hom_count > 0:
                    hh_ratio = round(float(het_count) / hom_count, 2)
                else:
                    hh_ratio = "NA"
            except ValueError:
                het_count = "NAO"
                hom_count = "NAO"
                hh_ratio = "NAO"
        else:
            het_count = "NAS"
            hom_count = "NAS"
            hh_ratio = "NAS"
            
        cnv_anno[key]['Het_count'] = het_count
        cnv_anno[key]['Hom_count'] = hom_count
        cnv_anno[key]['Het/Hom_ratio'] = hh_ratio
    
    '''
    vcf_reader1 = vcf.Reader(open(vcf_file, 'r'))
    cnv = pybedtools.BedTool(cnv_file)
    het_count = 0
    hom_count = 0
    
    for line in cnv:
        line[0] = line[0][3:]
        #bed_strng.append("\t".join(line[0:3]))
        key = ":".join(line[0:3])
        het_count = 0
        hom_count = 0
        for record in vcf_reader1.fetch(line[0], int(line[1]), int(line[2])):
            if record.num_het != 0 and len(record.FILTER) == 0:
            #if record.num_het != 0:
                het_count = het_count +1
            else:
                hom_count += 1
        if hom_count > 0:
            hh_ratio = round(float(het_count) / hom_count, 2)
        else:
            hh_ratio = 0
        
        cnv_res[key] = {}
        cnv_res[key]['CNV_st'] = line[3]
        cnv_res[key]['Het_count'] = het_count
        cnv_res[key]['Hom_count'] = hom_count
        cnv_res[key]['Het/Hom ratio'] = hh_ratio
        print "Region: %r CNV_STATE: %r Het_count: %r, Home_count: %r, Het/Hom ratio: %r" %(line[0:3], cnv_res[key]['CNV_st'], cnv_res[key]['Het_count'], cnv_res[key]['Hom_count'], cnv_res[key]['Het/Hom ratio'])
        '''

    return cnv_anno


def db_search(db_file, db_id, cnv_anno):
    for key in cnv_anno.keys():
        hit_count = {}
        hit_scores = []
        #db_id_count_lbl = db_id + "_count"
        #db_id_minmaxmedian_lbl = db_id + "_min_max_median"
        line = key.split(":")
        #print "%r" %line
        for row in db_file.fetch(str(line[0]), int(line[1]), int(line[2])):
            #print "%r" %row.split("\t")
            hit_count[row.split("\t")[5]] = 1
            hit_scores.append(float(row.split("\t")[4]))
        if len(hit_count) > 1:
            cnv_anno[key]['inDB_count'] = len(hit_count) - 1
            cnv_anno[key]['inDB_minmaxmedian'] = "|".join([str(min(hit_scores)), str(max(hit_scores)), str(np.median(hit_scores))])
            #print "hit_count: %r hit_scores: %r hit_scores (min): %r hit_scores (max): %r hit_score (median): %r" %(hit_count, hit_scores, min(hit_scores), max(hit_scores), np.median(hit_scores))
        else:
            cnv_anno[key]['inDB_count'] = "NA"
            cnv_anno[key]['inDB_minmaxmedian'] = "NA" 
            #print "no db_file hits"
    return cnv_anno


def array_search_run(cnv_file, db_file, freq = 0):
    db_filt_res = {}
    a_cnv = pybedtools.BedTool(cnv_file)
    b_db = pybedtools.BedTool(db_file)
    a_and_b = a_cnv.intersect(b_db, u=True)
    if freq == 0:
        for line in a_and_b:
            db_filt_res[":".join(line[0:3])] = "1"
    elif freq == 1:
        for line in a_and_b:
            if db_filt_res.get(":".join(line[0:3])):
                if db_filt_res[":".join(line[0:3])] < line[-1]:
                    db_filt_res[":".join(line[0:3])] = line[-1]
            else:
                db_filt_res[":".join(line[0:3])] = line[-1]
    return db_filt_res

def array_search(db_id, cnv_file, array_file, cnv_anno):
    db_filt_res = {}
    #db_filt_res = db_filt_run(cnv_file, db_file)
    a_cnv = pybedtools.BedTool(cnv_file)
    b_db = pybedtools.BedTool(array_file)
    a_and_b = a_cnv.intersect(b_db, u=True)

    for line in a_and_b:
        db_filt_res[":".join(line[0:3])] = "1"
    for key, v1 in db_filt_res.items():
        cnv_anno[key[3:]][db_id] = "PASS"

    return cnv_anno

'''
cnv_res_file = "/home/saneth/Documents/cnvFilt_proj/testSamples/A2013_exCopyDepth.tab"
vcf_file = pysam.TabixFile("/home/saneth/Documents/cnvFilt_proj/cnvScAn/PIDD/VCF_PIDD_txts_extracted.reformatted/P13-020_extracted.reformatted.tab.gz")
#vcf_file = "/home/saneth/Documents/cnvFilt_proj/testSamples/Berge-excap-A2013-Av5_all.filter.hgmd.vcf.gz"
db_file = pysam.TabixFile("/home/saneth/Documents/cnvFilt_proj/CNV_validation_test/excopydepth/inDB/ExCopyDepth_intersectBedC.tab.bed.gz")
db_id = "TTT"



cnv_anno = {}
cnv_anno, cnvs_ordered = read_cnvRes(cnv_res_file, cnv_anno)
cnv_anno = het_hom_count_2(cnv_anno, vcf_file)
db_search(db_file, "AA", cnv_anno)

for k in cnvs_ordered:
    print "key: %r Vale: %r" %(k,cnv_anno[k])
'''

'''
cnv_file = "../res-cnv/hChr"
vcf_file = "../indexed/all.filter.vcf.gz"
exonDef_file = "../resources/uniq.exon.havana_or_ensembl_gencode.v21.annotation.bed"
array_file = "../res-cnv/arr_h"
indb_file = "../res-cnv/indb_h"

cnv_res = {}

cnv_res = het_home_count(cnv_res, vcf_file, cnv_file)


#######
e_phastCon = pysam.TabixFile("../resources/PhastCon/phastConsElements100wayFormatted.bed.gz")
#######

b_gencode = pybedtools.BedTool("../res-cnv/1_gen_head.gtf")
c_cosCNV = pybedtools.BedTool("../resources/cosmicCompleteCNVs.tab")
d_conradCNV = pybedtools.BedTool("../resources/conrad.et.al.2010_Validated_CNVEs_v5_4Release.tab")
e_dgvCNV = pybedtools.BedTool("../resources/dgv_GRCh37_hg19_variants_2014-10-16.tab")
f_array = pybedtools.BedTool("../res-cnv/arr_h")
g_indb = pybedtools.BedTool("../res-cnv/indb_h")



cnv_res = filt_cnv(cnv_res, vcf_file, cnv_file, exonDef_file)
cnv_res= db_filt("aCGH", cnv_file, array_file, cnv_res)
cnv_res= db_filt("inDB", cnv_file, indb_file, cnv_res, 1)

cnv_res = phastCon_annotate(e_phastCon, cnv_res)


for k, v1 in cnv_res.items():
    line = k.split(":")
    lst_name = []
    exon_c = []    
    if cnv_res[k].get('Het_count'):
        line.append(str(cnv_res[k]['Het_count']))
    else:
        line.append("NA")
    if cnv_res[k].get('PVal'):
        line.append(str(cnv_res[k]['PVal']))
    else:
        line.append("NA")
    if cnv_res[k].get('Hom_count'):
        line.append(str(cnv_res[k]['Hom_count']))
    else:
        line.append("NA")
    if cnv_res[k].get('aCGH'):
        line.append(cnv_res[k]['aCGH'])
    else:
        line.append("NA")
    if cnv_res[k].get('inDB'):
        line.append(str(cnv_res[k]['inDB']))
    else:
        line.append("NA")
    # line.append("phastCon_count")
    line.append(str(cnv_res[k]['phastCon_count']))
    print "\t".join(line), len(line)
    pass
    
def filt_cnv(cnv_res, vcf_file, cnv_file, exonDef_file): #not used anymore
    exon_vcf = {}
    exon_cnv = {}
    cnv_vcf = {}
    plot_data = []
    bed_strng = []
    het_count = 0
    hom_count = 0
    
    vcf_reader1 = vcf.Reader(open(vcf_file, 'r'))
    a_exons = pybedtools.BedTool(exonDef_file)
    b_cnv = pybedtools.BedTool(cnv_file)
    
    for line in b_cnv:
        line[0] = line[0][3:]
        bed_strng.append("\t".join(line[0:3]))
        key = ":".join(line[0:3])
        for record in vcf_reader1.fetch(line[0], int(line[1]), int(line[2])):
            if record.num_het != 0:
                het_count = het_count +1
            else:
                hom_count += 1
        #print "Region: %r Het_count: %r, Home_count: %r" %(line[0:3], het_count, hom_count)
        cnv_res[key] = {}
        cnv_res[key]['Het_count'] = het_count
        cnv_res[key]['Hom_count'] = hom_count
    
    b_cnv_fmted = pybedtools.BedTool("\n".join(bed_strng), from_string=True)
    a_and_b = a_exons.intersect(b_cnv_fmted, wa=True, wb=True)
    #print a_and_b

    for line in a_and_b:
        key  = ":".join(line[0:3])
        val = ":".join(line[4:8])
        if exon_cnv.get(key):

            exon_cnv[key] = val
        else:
            exon_cnv[key] = val
    
    with open(exonDef_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            splt = line.split("\t")
            key  = ":".join(splt[0:3])
            for record in vcf_reader1.fetch(splt[0], int(splt[1]), int(splt[2])):
                #if record.num_het != 0 and len(record.FILTER) == 0:
                if record.num_het != 0:
                    if exon_vcf.get(key):
                        exon_vcf[key] = exon_vcf[key] + 1
                    else:
                        exon_vcf[key] = 1

    for k1, v1 in exon_vcf.items():
        plot_data.append(v1)

    for k1, v1 in exon_vcf.items():
        #print k1 + "**--**", v1
        if exon_cnv.get(k1):
            if cnv_vcf.get(exon_cnv[k1]):
                #print  k1, "******", exon_cnv[k1]
                t_arr = cnv_vcf[exon_cnv[k1]]
                t_arr.append(v1)
                cnv_vcf[exon_cnv[k1]] = t_arr
                pass
            else:
                #print  k1, "******", exon_cnv[k1]
                cnv_vcf[exon_cnv[k1]] = [v1] 
            pass
        pass

    for k1, v1 in cnv_vcf.items():
        #print k1 + "**^^^**", v1
        if len(v1) > 1:
            #print "****", mannwhitneyu(plot_data, v1, use_continuity=True)[1]
            cnv_res[k1]['PVal'] = mannwhitneyu(plot_data, v1, use_continuity=True)[1]
        #plot_data.append(v1)
        else:
            #print "---------------------------"
            cnv_res[k1]['PVal'] =  "NA"
        pass
    
    return(cnv_res)

'''
