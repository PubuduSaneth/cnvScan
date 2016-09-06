#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for f in $(find . -name "*.cnvScan.bed" -type f)
do
    wc -l $f
    python ${DIR}/../src/cnvScan_run.py $f ${f%.*}.out.tab ${DIR}/../data/DB/DB.txt.sorted.bed.gz
done
