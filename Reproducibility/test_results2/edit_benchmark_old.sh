#!/bin/bash


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = $3 - 1000
        new_end = $3
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/NK_B/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/NK_B/positive_data_filtered_right_padded.bed


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = center - 500
        new_end = center + 500
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/NK_B/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/NK_B/positive_data_filtered_center_padded.bed


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = $2
        new_end = $2 + 1000
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/NK_B/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/NK_B/positive_data_filtered_left_padded.bed


TARGET_COL=9

awk -v col="$TARGET_COL" -F'\t' '{
    line1 = $0
    val1 = $col
    getline line2 < "dataset/positive/NK_B/positive_data_filtered_right_padded.bed"
    split(line2, arr2, "\t")
    val2 = arr2[col]
    getline line3 < "dataset/positive/NK_B/positive_data_filtered_center_padded.bed"
    split(line3, arr3, "\t")
    val3 = arr3[col]
    max_val = val1 + 0
    best_line = line1
    if ((val2 + 0) > max_val) {
        max_val = val2 + 0
        best_line = line2
    }
    if ((val3 + 0) > max_val) {
        max_val = val3 + 0
        best_line = line3
    }
    print best_line
}' dataset/positive/NK_B/positive_data_filtered_left_padded.bed | awk -F'\t' -v OFS='\t' '{NF--; print}' > dataset/positive/NK_B/positive_data_filtered_edited.bed



awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = $3 - 1000
        new_end = $3
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/Monocyte_NK/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/Monocyte_NK/positive_data_filtered_right_padded.bed


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = center - 500
        new_end = center + 500
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/Monocyte_NK/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/Monocyte_NK/positive_data_filtered_center_padded.bed


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = $2
        new_end = $2 + 1000
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/Monocyte_NK/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/Monocyte_NK/positive_data_filtered_left_padded.bed




TARGET_COL=9

awk -v col="$TARGET_COL" -F'\t' '{
    line1 = $0
    val1 = $col
    getline line2 < "dataset/positive/Monocyte_NK/positive_data_filtered_right_padded.bed"
    split(line2, arr2, "\t")
    val2 = arr2[col]
    getline line3 < "dataset/positive/Monocyte_NK/positive_data_filtered_center_padded.bed"
    split(line3, arr3, "\t")
    val3 = arr3[col]
    max_val = val1 + 0  
    best_line = line1
    if ((val2 + 0) > max_val) { 
        max_val = val2 + 0
        best_line = line2 
    }
    if ((val3 + 0) > max_val) { 
        max_val = val3 + 0
        best_line = line3 
    }
    print best_line
}' dataset/positive/Monocyte_NK/positive_data_filtered_left_padded.bed | awk -F'\t' -v OFS='\t' '{NF--; print}' > dataset/positive/Monocyte_NK/positive_data_filtered_edited.bed


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = $3 - 1000
        new_end = $3
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/Monocyte_B/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/Monocyte_B/positive_data_filtered_right_padded.bed


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = center - 500
        new_end = center + 500
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/Monocyte_B/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/Monocyte_B/positive_data_filtered_center_padded.bed


awk 'BEGIN {OFS="\t"} {
    size = $3 - $2
    if (size < 1000) {
        center = int(($2 + $3) / 2)
        new_start = $2
        new_end = $2 + 1000
        if (new_start < 0) {
            new_start = 0
        }
        $2 = new_start
        $3 = new_end
    }
    print $0
}' dataset/positive/Monocyte_B/positive_data_filtered.bed | bedtools intersect -a - -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c > dataset/positive/Monocyte_B/positive_data_filtered_left_padded.bed




TARGET_COL=9

awk -v col="$TARGET_COL" -F'\t' '{
    line1 = $0
    val1 = $col
    getline line2 < "dataset/positive/Monocyte_B/positive_data_filtered_right_padded.bed"
    split(line2, arr2, "\t")
    val2 = arr2[col]
    getline line3 < "dataset/positive/Monocyte_B/positive_data_filtered_center_padded.bed"
    split(line3, arr3, "\t")
    val3 = arr3[col]
    max_val = val1 + 0  
    best_line = line1
    if ((val2 + 0) > max_val) { 
        max_val = val2 + 0
        best_line = line2 
    }
    if ((val3 + 0) > max_val) { 
        max_val = val3 + 0
        best_line = line3 
    }
    print best_line
}' dataset/positive/Monocyte_B/positive_data_filtered_left_padded.bed | awk -F'\t' -v OFS='\t' '{NF--; print}' > dataset/positive/Monocyte_B/positive_data_filtered_edited.bed
