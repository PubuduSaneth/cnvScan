#!/bin/bash

set -euf -o pipefail

mkdir -p DB
rm -f DB/DB.txt

i=0
for f in $(find . -name "*.cnvScan.bed" -type f)
do
    echo $f
    awk 'BEGIN{OFS="\t"}{printf "%s\tsample%03i\n", $0, $i}' $f >> DB/DB.txt
    i=$(($i+1))
done
sort -k1,1 -k2n,2 -k3n,3 DB/DB.txt > DB/DB.txt.sorted.bed
bgzip -f DB/DB.txt.sorted.bed
tabix -p bed DB/DB.txt.sorted.bed.gz
