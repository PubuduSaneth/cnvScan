#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for f in $(find . -name "*.cnvScan.bed" -type f)
do
    wc -l $f
    python ${DIR}/../src/cnvScan_run.py \
        -i $f \
        -o ${f%.*}.out.tab \
        -db ${DIR}/../data/DB/DB.txt.sorted.bed.gz \
        --resources ${DIR}/../resources/                    
done
