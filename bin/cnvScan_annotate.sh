#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python ${DIR}/../src/cnvScan_run.py \
    -i $@ \
    -o ${@%.*}.annotated.tsv \
    -db ${DIR}/../data/DB/DB.txt.sorted.bed.gz \
    --resources ${DIR}/../resources/
