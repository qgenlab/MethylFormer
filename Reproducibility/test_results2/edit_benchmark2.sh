#!/bin/bash

DATASETS=("NK_B" "Monocyte_NK" "Monocyte_B")
CPG_FILE="/mnt/analysis/derbelh/CpG_Without_Strand.bed"
TARGET_COL=9
MIN_CPG=5

for DS in "${DATASETS[@]}"; do
    DIR="dataset/positive/$DS"
    INPUT_BED="$DIR/positive_data_filtered.bed"
    awk 'BEGIN{OFS="\t"} {s=$3-$2; if(s<1000){start=$3-1000; end=$3; if(start<0)start=0; $2=start; $3=end} print $0}' "$INPUT_BED" | bedtools intersect -a - -b "$CPG_FILE" -c > "$DIR/right_tmp.bed"
    awk 'BEGIN{OFS="\t"} {s=$3-$2; if(s<1000){c=int(($2+$3)/2); start=c-500; end=c+500; if(start<0)start=0; $2=start; $3=end} print $0}' "$INPUT_BED" | bedtools intersect -a - -b "$CPG_FILE" -c > "$DIR/center_tmp.bed"
    awk 'BEGIN{OFS="\t"} {s=$3-$2; if(s<1000){start=$2; end=$2+1000; if(start<0)start=0; $2=start; $3=end} print $0}' "$INPUT_BED" | bedtools intersect -a - -b "$CPG_FILE" -c > "$DIR/left_tmp.bed"
    awk -v col="$TARGET_COL" -v n="$MIN_CPG" \
        -v orig_file="$INPUT_BED" \
        -v right_file="$DIR/right_tmp.bed" \
        -v center_file="$DIR/center_tmp.bed" \
        -F'\t' 'BEGIN {OFS="\t"} {
        getline line_orig < orig_file; split(line_orig, a_orig, "\t")
        orig_s = a_orig[2]; orig_e = a_orig[3]
        line1 = $0; val1 = $col
        getline line2 < right_file; split(line2, a2, "\t"); val2 = a2[col]
        getline line3 < center_file; split(line3, a3, "\t"); val3 = a3[col]
        max_val = val1 + 0; best_line = line1
        if ((val2 + 0) > max_val) { max_val = val2 + 0; best_line = line2 }
        if ((val3 + 0) > max_val) { max_val = val3 + 0; best_line = line3 }
        if (max_val >= n) {
            split(best_line, final_arr, "\t")
            final_arr[2] = orig_s
            final_arr[3] = orig_e
            for (i=1; i<col; i++) printf "%s%s", final_arr[i], (i==col-1 ? "" : OFS)
            print ""
        }
    }' "$DIR/left_tmp.bed" > "$DIR/positive_data_filtered_edited.bed"
    rm "$DIR/right_tmp.bed" "$DIR/center_tmp.bed" "$DIR/left_tmp.bed"
done
