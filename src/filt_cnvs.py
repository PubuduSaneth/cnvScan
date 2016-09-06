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


def db_search(db_file, cnv_anno):
    for key in cnv_anno.keys():
        hit_count = {}
        hit_scores = []
        line = key.split(":")
        for row in db_file.fetch(str(line[0]), int(line[1]), int(line[2])):
            hit_count[row.split("\t")[5]] = 1
            hit_scores.append(float(row.split("\t")[4]))
        if len(hit_count) > 1:
            cnv_anno[key]['inDB_count'] = len(hit_count) - 1
            cnv_anno[key]['inDB_minmaxmedian'] = "|".join([str(min(hit_scores)), str(max(hit_scores)), str(np.median(hit_scores))])
        else:
            cnv_anno[key]['inDB_count'] = "NA"
            cnv_anno[key]['inDB_minmaxmedian'] = "NA" 
    return cnv_anno


