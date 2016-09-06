#!/bin/bash

set -euf -o pipefail

for f in $(find . -name "*.cnvScan.bed" -type f)
do
    awk 'BEGIN{OFS="\t"}{print $0, substr(FILENAME,3)}' $f >> DB/DB.txt
done
sort -k1,1 -k2n,2 -k3n,3 DB/DB.txt > DB/DB.txt.sorted.bed
bgzip DB/DB.txt.sorted.bed
tabix -p bed DB/DB.txt.sorted.bed.gz
