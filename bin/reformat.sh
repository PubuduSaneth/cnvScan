#!/bin/bash

set -euf -o pipefail

for f in $(find . -name "*.cnv.bed")
do
    echo $f
    awk 'BEGIN \
    {
        OFS = "\t"
        cnv = ""
    }
    {
        if ($1 != "start.p") {
            if ($3 == "deletion") {cnv = 1}
            else if ($3 == "duplication") {cnv = 3}
            split($8, a, ":")
            split(a[2], b, "-")
            print substr(a[1], 4), b[1], b[2], cnv, $9
        }
    }' $f > ${f:2:7}'_cnvScan.ready.bed'
done
