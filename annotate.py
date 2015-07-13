#!/usr/bin/env python

import pybedtools
import pysam
import vcf
import re
import vcf
import numpy as np
from scipy.stats import mannwhitneyu
import filt_cnvs

def create_bedTools(cnv_file):

    cnv_bed = ""

    with open(cnv_file, 'r') as f:
        #next(f)
        for line in f:
            line = line.replace("\n","")
            line = line.replace("\r","")
            line = "chr"+line+'\n'
            cnv_bed = cnv_bed + line
    cnv = pybedtools.BedTool(cnv_bed,from_string=True)
    
    return cnv


def gencode_annotate(a_cnv, b_gencode, cnv_anno):
    
    utr_dict= {}
    transcript_dict= {}
    cnv_anno_list = []


    for line, feature in enumerate(a_cnv):
        f_id = ":".join(feature[0:3]) #f_id = ":".join(feature[0:3])[3:]
        cnv_anno_list.append(f_id)
    
    a_and_b = a_cnv.intersect(b_gencode, wa=True, wb=True)
    
    for line, feature in enumerate(a_and_b):
        feature_id = ":".join(feature[0:3])[3:]
        if cnv_anno.get(feature_id):
            pass
        else:
            cnv_anno[feature_id] = {}
        feature[-7] = feature[-7].replace('"','')
        if feature[-7] == 'gene':
            #print "111", feature
            g_cov = ""
            g_name = filter(lambda x: 'gene_name' in x, feature[-1].split(";"))[0].split(" ")[2].replace("\"","")
            g_type = filter(lambda x: 'gene_type' in x, feature[-1].split(";"))[0].split(" ")[2].replace("\"","")
            g_id = filter(lambda x: 'gene_id' in x, feature[-1].split(";"))[0].split(" ")[1].replace("\"","")
            
            if feature[-6] >= feature[1] and feature[-5] <= feature[2]:
                g_cov = "F"
            else:
                g_cov = "P"
            if cnv_anno[feature_id].get('gene_name'):
                pass
            else:
                cnv_anno[feature_id]['gene_name'] = {}
                cnv_anno[feature_id]['gene_type'] = {}
                cnv_anno[feature_id]['gene_id'] = {}
            
            cnv_anno[feature_id]['gene_name'][g_name] = g_cov
            cnv_anno[feature_id]['gene_type'][g_type] = 1
            cnv_anno[feature_id]['gene_id'][g_id] = 1
            
        elif feature[-7] == 'UTR':
            utr_gene = filter(lambda x: 'gene_name' in x, feature[-1].split(";"))[0].split(" ")[2].replace("\"","")
            if utr_dict.get(feature_id):
                pass
            else:
                utr_dict[feature_id] = {}
                
            if utr_dict[feature_id].get(utr_gene):
                utr_dict[feature_id][utr_gene] = ";".join([ utr_dict[feature_id][utr_gene], "-".join([feature[0], feature[-6], feature[-5]]) ])
            else:
                utr_dict[feature_id][utr_gene] = "-".join([feature[0], feature[-6], feature[-5]])
                
        elif feature[-7] == 'exon':
            feature[-1] = feature[-1].replace('"','')
            if cnv_anno[feature_id].get('exon'):
                pass
            else:
                cnv_anno[feature_id]['exon'] = {}
            if cnv_anno[feature_id]['exon'].get(filter(lambda x: 'transcript_id' in x, feature[-1].split(";"))[0].split(" ")[2]):
                cnv_anno[feature_id]['exon'][filter(lambda x: 'transcript_id' in x, feature[-1].split(";"))[0].split(" ")[2]] = cnv_anno[feature_id]['exon'][filter(lambda x: 'transcript_id' in x, feature[-1].split(";"))[0].split(" ")[2]] + 1
            else:
                cnv_anno[feature_id]['exon'][filter(lambda x: 'transcript_id' in x, feature[-1].split(";"))[0].split(" ")[2]] = 1
        
        elif feature[-7] == 'transcript':
            feature[-1] = feature[-1].replace('"','')
            #feature[-7] = feature[-7].replace('"','')
            transcript_gene = filter(lambda x: 'gene_name' in x, feature[-1].split(";"))[0].split(" ")[2].replace("\"","")
            if transcript_dict.get(feature_id):
                pass
            else:
                transcript_dict[feature_id] = {}

            if transcript_dict[feature_id].get(transcript_gene):
                transcript_dict[feature_id][transcript_gene] = ";".join([ transcript_dict[feature_id][transcript_gene], filter(lambda x: 'transcript_id' in x, feature[-1].split(";"))[0].split(" ")[2] ])
            else:
                transcript_dict[feature_id][transcript_gene] = {}
                transcript_dict[feature_id][transcript_gene] = filter(lambda x: 'transcript_id' in x, feature[-1].split(";"))[0].split(" ")[2] #"-".join([feature[0], feature[-6], feature[-5]])
    
    for k in cnv_anno.keys():
        if utr_dict.get(k):
            for g in utr_dict[k].keys():
                if cnv_anno[k]['gene_name'].get(g):
                    if cnv_anno[k]['gene_name'][g] == "P":
                        if cnv_anno[k].get('UTR'):
                            cnv_anno[k]['UTR'] = "|".join([cnv_anno[k]['UTR'], utr_dict[k][g]])
                        else:
                            cnv_anno[k]['UTR'] = utr_dict[k][g]
        '''if cnv_anno[k].get('UTR'):
            pass
        else:
            cnv_anno[k]['UTR'] = "NA"
        #print "UUUUUUUUUUUUUTTTTTTTTTTTTRRRRRRRRRRRR", k, cnv_anno[k]['UTR']'''
    

    for k in cnv_anno.keys():
        if transcript_dict.get(k):
            for g in transcript_dict[k].keys():
                if cnv_anno[k]['gene_name'].get(g):
                    if cnv_anno[k]['gene_name'][g] == "P":
                        if cnv_anno[k].get('transcript'):
                            cnv_anno[k]['transcript'] = "|".join([cnv_anno[k]['transcript'], transcript_dict[k][g]])
                        else:
                            cnv_anno[k]['transcript'] = transcript_dict[k][g]
                        if cnv_anno[k].get('exon'):
                            if cnv_anno[k].get('exon_count'):
                                if cnv_anno[k]['exon'].get(trans_id):
                                    cnv_anno[k]['exon_count'][trans_id] = cnv_anno[k]['exon'][trans_id]
                            else:
                                cnv_anno[k]['exon_count'] = {}
                                for trans_id in transcript_dict[k][g].split(";"):
                                    #print trans_id
                                    if cnv_anno[k]['exon'].get(trans_id):
                                        cnv_anno[k]['exon_count'][trans_id] = cnv_anno[k]['exon'][trans_id]

        '''if cnv_anno[k].get('transcript'):
            pass
        else:
            cnv_anno[k]['transcript'] = "NA"
        print "transcripttranscripttranscript",k, cnv_anno[k]['transcript']'''

    '''for k in cnv_anno.keys():
        if cnv_anno[k].get('exon'):
            cnv_anno[k]['exon_count'] = 0
            for trans_id in cnv_anno[k]['exon'].keys():
                if cnv_anno[k]['exon'][trans_id] > cnv_anno[k]['exon_count']:
                    cnv_anno[k]['exon_count'] = cnv_anno[k]['exon'][trans_id]
        if cnv_anno[k].get('exon_count'):
            pass
        else:
            cnv_anno[k]['exon_count'] = "NA"
        print "exon_countexon_countexon_count",k, cnv_anno[k]['exon_count'] '''

    return cnv_anno

def sanger_annotate(a_cnv, c_conradCNV, cnv_anno):
    a_and_c = a_cnv.intersect(c_conradCNV, wa=True, wb=True)

    for line, feature in enumerate(a_and_c):
        feature_id = ":".join(feature[0:3])[3:]
        if cnv_anno.get(feature_id):
            pass
        else:
            cnv_anno[feature_id] = {}    
        if cnv_anno[feature_id].get('Sanger_HiRes_CNV'):
            cnv_anno[feature_id]['Sanger_HiRes_CNV'] = cnv_anno[feature_id]['Sanger_HiRes_CNV'] + 1
        else:
            cnv_anno[feature_id]['Sanger_HiRes_CNV'] = 1
    return cnv_anno

def dgv_annotate(a_cnv, d_dgvCNV, cnv_anno):
    a_and_e = a_cnv.intersect(d_dgvCNV, wa=True, wb=True)   

    for line, feature in enumerate(a_and_e):
        #print line, feature[0:8]
        feature_id = ":".join(feature[0:3])[3:]
        if cnv_anno.get(feature_id):
            pass
        else:
            cnv_anno[feature_id] = {}    
        if cnv_anno[feature_id].get('DGV_CNV'):
            cnv_anno[feature_id]['DGV_CNV'] = cnv_anno[feature_id]['DGV_CNV'] + 1
        else:
            cnv_anno[feature_id]['DGV_CNV'] = 1
        if cnv_anno[feature_id].get('DGV_VarType'):
            if filter(lambda x: feature[-3] in x, cnv_anno[feature_id]['DGV_VarType'].split(";")):
                pass
            else:
                cnv_anno[feature_id]['DGV_VarType'] = ";".join([ cnv_anno[feature_id]['DGV_VarType'], feature[-3] ])
        else:
            cnv_anno[feature_id]['DGV_VarType'] = feature[-3]
        if cnv_anno[feature_id].get('DGV_VarSubType'):
            if filter(lambda x: feature[-2] in x, cnv_anno[feature_id]['DGV_VarSubType'].split(";")):
                pass
            else:
                cnv_anno[feature_id]['DGV_VarSubType'] = ";".join([ cnv_anno[feature_id]['DGV_VarSubType'], feature[-2] ])
        else:
            cnv_anno[feature_id]['DGV_VarSubType'] = feature[-2]
        if cnv_anno[feature_id].get('DGV_PUBMEDID'):
            if filter(lambda x: feature[-1] in x, cnv_anno[feature_id]['DGV_PUBMEDID'].split(";")):
                pass
            else:
                cnv_anno[feature_id]['DGV_PUBMEDID'] = ";".join([ cnv_anno[feature_id]['DGV_PUBMEDID'], feature[-1] ])
        else:
            cnv_anno[feature_id]['DGV_PUBMEDID'] = feature[-1]
    return cnv_anno

def dgvFilt_annotate(d_dgvFiltsCNV_l, cnv_anno, level):
    for k in cnv_anno.keys():
        dgvFilt_count = 0
        line = k.split(":")
        dgv_popFreq = []
        lbl_count = level + "_count"
        lbl_pop = level + "_popFreq"
        
        #chrom = 'chr' + line[0]
        try:
            for row in d_dgvFiltsCNV_l.fetch(line[0], int(line[1]), int(line[2])):
                dgvFilt_count = dgvFilt_count  + 1
                dgv_popFreq.append(row.split("\t")[4])
            if len(dgv_popFreq) > 0:
                cnv_anno[k][lbl_count] = dgvFilt_count
                cnv_anno[k][lbl_pop] = "|".join(dgv_popFreq)
                #print "11111111||| k: %s >> %s: %d | %s: %s" %(k, lbl_count, cnv_anno[k][lbl_count], lbl_pop, cnv_anno[k][lbl_pop])
            else:
                cnv_anno[k][lbl_count] = "NA"
                cnv_anno[k][lbl_pop] = "NA"
        except ValueError:
            pass
    return cnv_anno

def phastCon_annotate(e_phastCon, cnv_anno):
    for k in cnv_anno.keys():
        phastConEle_count = 0
        line = k.split(":")
        phastCon_lod = []
        chrom = 'chr' + line[0]
        #chrom = line[0]
        for row in e_phastCon.fetch(chrom, int(line[1]), int(line[2])):
            phastConEle_count = phastConEle_count  + 1
            phastCon_lod.append(int(row.split("\t")[3][4:]))
        if len(phastCon_lod) > 0:
            cnv_anno[k]['phastCon_count'] = phastConEle_count
            cnv_anno[k]['phastCon_min_max'] = ":".join([str(min(phastCon_lod)), str(max(phastCon_lod))])
        else:
            cnv_anno[k]['phastCon_count'] = "NA"
            cnv_anno[k]['phastCon_min_max'] = "NA"
    return cnv_anno
    pass

def geneticIntolarance_annotate(i_genIntol_file, cnv_anno):
    genticIntol_score = {}
    with open(i_genIntol_file, 'r') as f:
        for line in f:
            line = line.replace('\n', '')
            line = line.replace('\r', '')
            l = line.split("\t")
            genticIntol_score[l[0]] = l[1]    
    
    for k in cnv_anno.keys():
        #print "******", k
        if cnv_anno[k].get('gene_name'):
            for g in cnv_anno[k]['gene_name'].keys():
                if genticIntol_score.get(g):
                    #print "GGGG",g,k, genticIntol_score[g]
                    if cnv_anno[k].get('GenInTolScore'):
                        #print "000", cnv_anno[k]['GenInTolScore'], genticIntol_score[g]
                        #print cnv_anno[k]['GenInTolScore'], genticIntol_score[g]
                        cnv_anno[k]['GenInTolScore'] = "|".join([cnv_anno[k]['GenInTolScore'],genticIntol_score[g]])
                        #print "222: %r %r " %(k, cnv_anno[k]['GenInTolScore'])
                    else:
                        cnv_anno[k]['GenInTolScore'] = genticIntol_score[g]
                        #print "111: %r %r " %(k, cnv_anno[k]['GenInTolScore'])
                        #print "111", cnv_anno[k]['GenInTolScore'], genticIntol_score[g]
                        #print "111", cnv_anno[k]['GenInTolScore'], genticIntol_score[g]
                else:
                    if cnv_anno[k].get('GenInTolScore'):
                        pass
                    else:
                        cnv_anno[k]['GenInTolScore'] = "NA"
        else:
            cnv_anno[k]['GenInTolScore'] = "NA"
        #print "000000000000", cnv_anno[k]['GenInTolScore'] 
    return cnv_anno
    pass

def haploIdx_annotate(f_haploIdx, cnv_anno):
    for k in cnv_anno.keys():
        haploIdx_count = 0
        line = k.split(":")
        haploIdx_percentage = []
        haploIdx_score = []
        chrom = 'chr' + line[0]
        for row in f_haploIdx.fetch(chrom, int(line[1]), int(line[2])):
            haploIdx_count = haploIdx_count  + 1
            haploIdx_percentage.append(row.split("\t")[4][:-1])
            haploIdx_score.append(row.split("\t")[5])
        if len(haploIdx_percentage) > 0:
            cnv_anno[k]['haploIdx_count'] = haploIdx_count
            cnv_anno[k]['haploIdx_score'] = "|".join([":".join(haploIdx_percentage), ":".join(haploIdx_score)])
        else:
            cnv_anno[k]['haploIdx_count'] = "NA"
            cnv_anno[k]['haploIdx_score'] = "NA"
    return cnv_anno
    pass

def del1000g_annotate(g_del1000g_delFile, cnv_anno):
    for k in cnv_anno.keys():
        del_1000g_count = 0
        line = k.split(":")
        for row in g_del1000g_delFile.fetch(line[0][3:], int(line[1]), int(line[2])):
            del_1000g_count = del_1000g_count  + 1
        if del_1000g_count > 0:
            #print "Key: %r Count: %r Max: %r Min: %r" %(line, phastConEle_count, max(phastCon_lod), min(phastCon_lod))
            cnv_anno[k]['1000G_Del_count'] = del_1000g_count
        else:
            cnv_anno[k]['1000G_Del_count'] = "NA"
    return cnv_anno
    pass

def dup1000g_annotate(h_dup1000g_delFile, cnv_anno):
    for k in cnv_anno.keys():
        dup_1000g_count = 0
        line = k.split(":")
        if re.search(r'Y',line[0]):
            pass
        else:
            for row in h_dup1000g_delFile.fetch(line[0][3:], int(line[1]), int(line[2])):
                dup_1000g_count = dup_1000g_count  + 1
        if dup_1000g_count > 0:
            cnv_anno[k]['1000G_Dup_count'] = dup_1000g_count
        else:
            cnv_anno[k]['1000G_Dup_count'] = "NA"
    return cnv_anno
    pass

def clinVar_annotate(i_clinVar_reader, cnv_anno):

    for k in cnv_anno.keys():
        line = k.split(":")
        clindbn = {}
        clinhgvs = {}
        for record in i_clinVar_reader.fetch(line[0][3:], int(line[1]), int(line[2])):
            if len(record.REF) > 1 or len(record.ALT[0]) > 1:
                if re.search(r'4|5|6|7', record.INFO['CLNSIG'][0]):
                    #print "00000000000000000"
                    #print record,"\t",record.INFO['CLNHGVS'], record.INFO['CLNDSDB'], record.INFO['CLNSIG'], record.INFO['CLNDBN'], len(record.INFO['CLNDBN'])
                    for clin_disease in record.INFO['CLNDBN']:
                        for clin_disease_split in clin_disease.split("|"):
                            if clin_disease_split != 'not_provided':
                                clindbn[clin_disease_split] = 1
                            
                    for clinhgvs_name in record.INFO['CLNHGVS']:
                        clinhgvs[clinhgvs_name] = 1
        if len(clindbn.keys()) > 0:
            cnv_anno[k]['clindbn'] = "|".join(clindbn.keys())
        else:
            cnv_anno[k]['clindbn'] = "NA"
        if len(clinhgvs.keys()) > 0:
            cnv_anno[k]['clinhgvs'] = "|".join(clinhgvs.keys())
        else:
            cnv_anno[k]['clinhgvs'] = "NA"
    return cnv_anno

def omim_annotate(j_omim_file, cnv_anno):
    omim_morbidMap = {}
    with open(j_omim_file, 'r') as f:
        for line in f:
            line = line.replace('\n', '')
            l = line.split("\t")
            omim_morbidMap[l[0]] = l[1]

    for k in cnv_anno.keys():
        if cnv_anno[k].get('gene_name'):
            for g in cnv_anno[k]['gene_name'].keys():
                if omim_morbidMap.get(g):
                    if cnv_anno[k].get('OMIM'):
                        cnv_anno[k]['OMIM'] = "|".join([cnv_anno[k]['OMIM'],omim_morbidMap[g]])
                    else:
                        cnv_anno[k]['OMIM'] = omim_morbidMap[g]
                else:
                    if cnv_anno[k].get('OMIM'):
                        pass
                    else:
                        cnv_anno[k]['OMIM'] = "NA"
        else:
            cnv_anno[k]['OMIM'] = "NA"    
    return cnv_anno

def devDisorder_annotate(h_devDis_file, cnv_anno):
    devDisorder = {}
    with open(h_devDis_file, 'r') as f:
        for line in f:
            line = line.replace('\n', '')
            l = line.split("\t")
            devDisorder[l[0]] = l[1].split("|")
            #print l[0], devDisorder[l[0]]

    for k in cnv_anno.keys():
        if cnv_anno[k].get('gene_name'):
            cnv_anno_devDis = {}
            cnv_anno_devDis['devDis_mutConseq'] = {}
            cnv_anno_devDis['devDis_disName'] = {}
            cnv_anno_devDis['devDis_pubmedID'] = {}
            
            for g in cnv_anno[k]['gene_name'].keys():
                if devDisorder.get(g):
                    cnv_anno_devDis['devDis_mutConseq'][devDisorder[g][0]] = 1                    
                    cnv_anno_devDis['devDis_disName'][devDisorder[g][1]] = 1
                    cnv_anno_devDis['devDis_pubmedID'][devDisorder[g][2]] = 1
            
            if cnv_anno_devDis['devDis_mutConseq']:
                cnv_anno[k]['devDis_mutConseq'] = ";".join(cnv_anno_devDis['devDis_mutConseq'].keys())
            else:
                cnv_anno[k]['devDis_mutConseq'] =  "NA"
            if cnv_anno_devDis['devDis_disName']:
                cnv_anno[k]['devDis_disName'] = ";".join(cnv_anno_devDis['devDis_disName'].keys())
            else:
                cnv_anno[k]['devDis_disName'] = "NA"
            if cnv_anno_devDis['devDis_pubmedID']:
                cnv_anno[k]['devDis_pubmedID'] = ";".join(cnv_anno_devDis['devDis_pubmedID'].keys())
            else:
                cnv_anno[k]['devDis_pubmedID'] = "NA"

            #print k,cnv_anno[k]['gene_name'],cnv_anno[k]['devDis_mutConseq'], cnv_anno[k]['devDis_disName'], cnv_anno[k]['devDis_pubmedID']
            
        else:
            cnv_anno[k]['devDis_mutConseq'] =  "NA"
            cnv_anno[k]['devDis_disName'] =  "NA"
            cnv_anno[k]['devDis_pubmedID'] =  "NA"
            #print k,cnv_anno[k]['devDis_mutConseq'], cnv_anno[k]['devDis_disName'], cnv_anno[k]['devDis_pubmedID']
    
    return cnv_anno