#!/bin/bash


RUN_DIR=$PWD

bedtools intersect -a $RUN_DIR/results/NK_B_benchmark_3_tools.csv -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -wo | awk 'BEGIN {OFS="\t"} {dmr=$1"\t"$2"\t"$3; count[dmr]++} END {for (i in count) { split(i, coords, "\t"); start=coords[2]; end=coords[3]; size=(end - start + 1); ratio=count[i]/size; print i, count[i], size, ratio }}' | sort -k1,1 -k2,2n > $RUN_DIR/results/NK_B_benchmark_3_tools_nbr_CpG.bed

bedtools intersect -a $RUN_DIR/results/NK_Monocyte_benchmark_3_tools.csv -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -wo | awk 'BEGIN {OFS="\t"} {dmr=$1"\t"$2"\t"$3; count[dmr]++} END {for (i in count) { split(i, coords, "\t"); start=coords[2]; end=coords[3]; size=(end - start + 1); ratio=count[i]/size; print i, count[i], size, ratio }}' | sort -k1,1 -k2,2n > $RUN_DIR/results/NK_Monocyte_benchmark_3_tools_nbr_CpG.bed

bedtools intersect -a $RUN_DIR/results/B_Monocyte_benchmark_3_tools.csv -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -wo | awk 'BEGIN {OFS="\t"} {dmr=$1"\t"$2"\t"$3; count[dmr]++} END {for (i in count) { split(i, coords, "\t"); start=coords[2]; end=coords[3]; size=(end - start + 1); ratio=count[i]/size; print i, count[i], size, ratio }}' | sort -k1,1 -k2,2n > $RUN_DIR/results/B_Monocyte_benchmark_3_tools_nbr_CpG.bed


##############################################################################


awk 'BEGIN {OFS="\t"} {print $11, $12, $13}' $RUN_DIR/results/B_NK/dl_dmr_035_generate_DMR_0.005.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -wo | awk 'BEGIN {OFS="\t"} {dmr=$1"\t"$2"\t"$3; count[dmr]++} END {for (i in count) { split(i, coords, "\t"); start=coords[2]; end=coords[3]; size=(end - start + 1); ratio=count[i]/size; print i, count[i], size, ratio }}' | sort -k1,1 -k2,2n > $RUN_DIR/results/NK_B_dl_model_nbr_CpG.bed


awk 'BEGIN {OFS="\t"} {print $11, $12, $13}' $RUN_DIR/results/B_Monocytes/dl_dmr_035_generate_DMR_0.005.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -wo | awk 'BEGIN {OFS="\t"} {dmr=$1"\t"$2"\t"$3; count[dmr]++} END {for (i in count) { split(i, coords, "\t"); start=coords[2]; end=coords[3]; size=(end - start + 1); ratio=count[i]/size; print i, count[i], size, ratio }}' | sort -k1,1 -k2,2n > $RUN_DIR/results/B_Monocytes_dl_model_nbr_CpG.bed



awk 'BEGIN {OFS="\t"} {print $11, $12, $13}' $RUN_DIR/results/NK_Monocytes/dl_dmr_035_generate_DMR_0.005.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -wo | awk 'BEGIN {OFS="\t"} {dmr=$1"\t"$2"\t"$3; count[dmr]++} END {for (i in count) { split(i, coords, "\t"); start=coords[2]; end=coords[3]; size=(end - start + 1); ratio=count[i]/size; print i, count[i], size, ratio }}' | sort -k1,1 -k2,2n > $RUN_DIR/results/NK_Monocytes_dl_model_nbr_CpG.bed
