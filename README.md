# cnvScan
CNV screening and annotation tool

With advances in next generation sequencing technology and analysis methods, single nucleotide variants (SNVs) and indels can be detected with high sensitivity and specificity in exome sequencing data. While recent studies have demonstrated the ability to detect disease-causing copy number variants (CNVs), exonic CNV prediction programs have shown high false positive (FP) CNV counts. The high false discovery rate is the major limiting factor for the applicability of these programs in clinical studies. Thus, we developed a tool (cnvScan) to improve the clinical utility of computational CNV prediction by reducing the false positive count and providing annotations.

cnvScan consists of two stages: CNV screening and CNV annotation. The screening stage utilizes CNV quality scores reported from prediction programs and further evaluates the quality of CNV calls by using an in-house CNV database. The annotation stage is designed to annotate predicted CNVs with functionally and clinically significant information using different source datasets. In addition, known CNVs within predicted variants were annotated using publicly available CNV datasets. CNV quality scores and annotated information reported from cnvScan can be used to filter high-quality rare CNVs from computational CNV predictions.

For more information:
http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2374-2
