#!/bin/bash

BAM=$1

BAM_NAME=$( echo $BAM | xargs basename | awk -F. '{print $1}' | awk 'gsub("_RMLizer","",$1)' )

echo $BAM_NAME
echo "--> Running BedTools..."

bedtools bamtobed -i ${BAM_NAME}/${BAM_NAME}_RMLizer.bam -splitD  \
| gawk 'BEGIN{FS=OFS="\t"}; {gsub("/[12]$", "", $4); print $1, $2, $3, $4}' \
| sort -k4 \
| gawk 'BEGIN{FS=OFS="\t"}; {if (NR == 1){chr = $1; strt = $2; end = $3; rname = $4; mate = "A"} else if ($4 == rname){mate = "B"; if ($2 > end){print chr, strt, end, rname "_" mate; chr = $1; strt = $2; end = $3; } else {if ($3 > end){end = $3; }}} else if ($4 != rname){mate = "A"; print chr, strt, end, rname "_" mate; chr = $1; strt = $2; end = $3; rname = $4; }} END{print chr, strt, end, rname "_" mate}' \
| sort-bed --max-mem 8G --tmpdir ${BAM_NAME} - \
> ${BAM_NAME}/${BAM_NAME}_RMLizer_bamtobed_splitD_sortBed_merged.bed

echo "--> Half way through..."

cat ${BAM_NAME}/${BAM_NAME}_mm_pos_sorted.txt \
| sort-bed --max-mem 8G --tmpdir ${BAM_NAME} - \
| uniq \
| gawk 'BEGIN{FS=OFS="\t"}; {if ($8 == "F") {print $1, $2, $3, $4, ".", $5, $6, $7}}' \
| bedtools coverage -a stdin -b ${BAM_NAME}/${BAM_NAME}_RMLizer_bamtobed_splitD_sortBed_merged.bed -counts -sorted \
> ${BAM_NAME}/${BAM_NAME}_coverage.bed

echo "--> Done!"