#!/bin/bash

set -euf -o pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python ${DIR}/../src/cnvScan_VarFilt.py \
    -i $@ \
    -o ${@%.*}.filtered.tsv \
    --score 10 \
    --genelist genelist.tsv.gz
