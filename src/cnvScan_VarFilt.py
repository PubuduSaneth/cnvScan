#!/usr/bin/env python

import pysam
import sys
import argparse
from collections import OrderedDict

class cnv_filter(object):

    def __init__(self, input_file, output_file, score, genelist):
        self.input = input_file
        self.output = output_file
        self.score = score
        self.genelist = genelist
        self.process()

    def filter_line(self, row):
        keep = True
        keep = (float(row['default_score']) >= self.score)
        keep &= ((row['inDBScore_MinMaxMedian'] == 'NA') or ((row['inDBScore_MinMaxMedian'] != 'NA') and (float(row['inDBScore_MinMaxMedian'].split("|")[2]) >= self.score)))
        keep &= (row['DGV_Stringency2_count'] == "NA" or row['DGV_Stringency12_count'] == "NA")
        keep &= (row['1000g_del'] == "NA" or row['1000g_ins'] == "NA")
        if self.genelist:
            keep &= row.get('PIDD_GENE', 'NA') != 'NA'
        return keep

    def process(self):
        if self.genelist:
            gene_list = pysam.TabixFile(self.genelist)

        with open(self.input, 'r') as fin:
            with open(self.output, 'w') as fout:
                for line in fin:

                    line = line.strip()
                    if line.startswith('#chr'):
                        header = line.split('\t')
                        if self.genelist:
                            header += ['PIDD_GENE', 'Inheritance', 'Phenotype']
                        fout.write('\t'.join(header) + '\n')
                        continue
                    elif line.startswith('##'):
                        continue

                    try:
                        row = OrderedDict(zip(header, line.split('\t')))
                    except:
                        continue

                    if self.genelist:
                        PIDD_GENE = []
                        Inheritance = []
                        Phenotype = []
                        try:
                            for genepanel_line in gene_list.fetch(row['#chr'],
                                                                  int(row['start']),
                                                                  int(row['end']),
                                                                  parser=pysam.asTuple()):
                                PIDD_GENE += [genepanel_line[3]]
                                Inheritance += [genepanel_line[4]]
                                Phenotype += [genepanel_line[5]]
                            row['PIDD_GENE'] = "|".join(PIDD_GENE) if PIDD_GENE else 'NA'
                            row['Inheritance'] = "|".join(Inheritance) if Inheritance else 'NA'
                            row['Phenotype'] ="|".join(Phenotype) if Phenotype else 'NA'
                        except ValueError:
                            pass

                    if self.filter_line(row):
                        fout.write('\t'.join(row.values()) + '\n')


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
