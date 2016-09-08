#!/usr/bin/env python

import pysam
import sys
import argparse

'''
Run command
cnvVarFilt_v3.py quality_score path_to_cnvScan_files path_to_geneList

Field IDs
	0 . chr
	1 . start
	2 . end
	3 . cnv_state
	4 . score
	5 . len
	6 . inDB_count
	7 . inDB_MinMaxMedian
	8 . gene_name
	9 . gene_type
	10 . gene_id
	11 . exon_count
	12 . UTR
	13 . transcript
	14 . phastConElement_count
	15 . phastConElement_minMax
	16 . haplo_insufIdx_count
	17 . haplo_insufIdx_score
	18 . Gene_intolarance_score
	19 . sanger_cnv
	20 . dgv_cnv
	21 . dgv_varType
	22 . dgv_varSubType
	23 . dgv_pubmedId
	24 . DGV_Stringency2_count
	25 . DGV_Stringency2_PopFreq
	26 . DGV_Stringency12_count
	27 . DGV_Stringency12_popFreq
	28 . 1000g_del
	29 . 1000g_ins
	30 . omim_morbidMap
	31 . ddd_mutConsequence
	32 . ddd_diseaseName
	33 . ddd_pubmedId
	34 . clinVar_disease
	35 . hgvs_varName

'''

class cnv_filter(object):

    def __init__(self, input_file, output_file, score, genelist):
        self.input = input_file
        self.output = output_file
        self.score = score
        self.genelist = genelist
        self.process()


    def process(self):
		if self.genelist:
			gene_list = pysam.TabixFile(self.genelist)

		print """Filter protocol
		            4 score > %d and
		            7 inDB_MinMaxMedian.split("|")[2] > %d and
		            24 DGV_Stringency2_count == "NA" OR 26 DGV_Stringency12_count == "NA" and
		            28 1000g_del == "NA"
		            29 1000g_ins == "NA" """ %(self.score, self.score)

		out_file = open(self.output, 'w')

		out_file.write ("""Filtration protocol
		            Default score > %d and
		            CNVQ: inDB_MinMaxMedian.split("|")[2] > %d and
		            DGV CNVs: DGV_Stringency2_count == "NA" OR 26 DGV_Stringency12_count == "NA" and
		            1000 genome deletion: 1000g_del == "NA"
		            1000 genome insertion: 1000g_ins == "NA"\n""" %(self.score, self.score))

		if self.genelist:
			out_file.write("\t".join(['chr', 'start', 'end', 'cnv_state', 'default_score', 'len', 'inDB_count', 'inDBScore_MinMaxMedian', 'gene_name', 'gene_type', 'gene_id', 'exon_count', 'UTR', 'transcript', 'phastConElement_count', 'phastConElement_minMax', 'haplo_insufIdx_count', 'haplo_insufIdx_score', 'Gene_intolarance_score', 'sanger_cnv', 'dgv_cnv', 'dgv_varType', 'dgv_varSubType', 'dgv_pubmedId', 'DGV_Stringency2_count', 'DGV_Stringency2_PopFreq', 'DGV_Stringency12_count', 'DGV_Stringency12_popFreq', '1000g_del', '1000g_ins', 'omim_morbidMap', 'ddd_mutConsequence', 'ddd_diseaseName', 'ddd_pubmedId', 'clinVar_disease', 'hgvs_varName', 'PIDD_GENE', 'Inheritnce', 'Phenotype']) + "\n")
		else:
			out_file.write("\t".join(['chr', 'start', 'end', 'cnv_state', 'default_score', 'len', 'inDB_count', 'inDBScore_MinMaxMedian', 'gene_name', 'gene_type', 'gene_id', 'exon_count', 'UTR', 'transcript', 'phastConElement_count', 'phastConElement_minMax', 'haplo_insufIdx_count', 'haplo_insufIdx_score', 'Gene_intolarance_score', 'sanger_cnv', 'dgv_cnv', 'dgv_varType', 'dgv_varSubType', 'dgv_pubmedId', 'DGV_Stringency2_count', 'DGV_Stringency2_PopFreq', 'DGV_Stringency12_count', 'DGV_Stringency12_popFreq', '1000g_del', '1000g_ins', 'omim_morbidMap', 'ddd_mutConsequence', 'ddd_diseaseName', 'ddd_pubmedId', 'clinVar_disease', 'hgvs_varName']) + "\n")


		def filt_line(f_line):
			r_line = []

			if float(f_line[4]) > self.score and \
				(f_line[7] == "NA" or ( not f_line[7] == "NA" and float(f_line[7].split("|")[2]) > self.score)) and \
				(f_line[24] == "NA" or f_line[26] == "NA") and \
				(f_line[28] == "NA" or f_line[29] == "NA"):
				r_line = f_line
			return r_line

		with open(self.input, "r") as f1:

			next(f1)
			for line in f1:

				g_name = []
				inheritance = []
				phenotype = []
				line = line.replace("\n","")
				line = line.replace("\r","")
				line = line.split("\t")
				if self.genelist:
					try:
						for row in gene_list.fetch(str(line[0]), int(line[1]), int(line[2])):
							row = row.split("\t")
							g_name.append(row[3])
							inheritance.append(row[4])
							phenotype.append(row[5])
						if len(g_name) >0 :
							line.extend([ "|".join(g_name), "|".join(inheritance), "|".join(phenotype)])
							p_line = filt_line(line)
							if len(p_line) > 0:
								out_file.write("\t".join(p_line)+"\n")
					except ValueError:
						pass
				else:
					p_line = filt_line(line)
					if len(p_line) > 0:
						out_file.write("\t".join(p_line)+"\n")
		out_file.close()


# ============================================================
def main():
        '''
        cnvScan_VarFilt
        '''
        parser = argparse.ArgumentParser(
            description='Filter annotated CNV prediction based on rules')
        parser.add_argument('-i', '--input', type=str, required=True,
                            help='Input CNV bed file')
        parser.add_argument('-o', '--output', type=str, required=True,
                            help='Output CNV bed file')
        parser.add_argument('-s', '--score', type=float, required=False,
							default=10.0,
                            help='Quality score threshold')
        parser.add_argument('-g', '--genelist', type=str, required=False,
                            help='Gene list')
        args = parser.parse_args()

        filtered = cnv_filter(args.input, args.output, args.score, args.genelist)


# ============================================================
if __name__ == "__main__":
    sys.exit(main())
