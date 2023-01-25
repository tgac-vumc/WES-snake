#!/bin/bash

input=$1
sample=$2
outdir=$3
output=$4

samtools view -F 0x0404 "${input}" | awk "BEGIN {FS=\"\t\";uni=0;q1=0;q37=0} \
        {if (\$2 ~ /[^d]/ && \$5 >= 1) q1 += 1; \
        if (\$2 ~ /[^d]/ && \$5 >= 37) q37 += 1 } \
        END {print NR > \"${outdir}/${sample}.reads.unique\";\
        print q1 > \"${outdir}/${sample}.reads.q1\";\
        print q37 > \"${outdir}/${sample}.reads.q37\"}"

samtools idxstats "${input}" | awk "BEGIN {FS=\"\t\";aligned=0;unaligned=0} {aligned += \$3; unaligned += \$4} END {print aligned > \"${outdir}/${sample}.reads.aligned\"; print aligned + unaligned > \"${output}\" } "
