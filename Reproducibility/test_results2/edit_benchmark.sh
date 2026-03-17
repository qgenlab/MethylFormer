#!/bin/bash

# 1. Define your datasets and your CpG background file
DATASETS=("NK_B" "Monocyte_NK" "Monocyte_B")
CPG_FILE="/mnt/analysis/derbelh/CpG_Without_Strand.bed"
TARGET_COL=9

# 2. Loop over each dataset
for DS in "${DATASETS[@]}"; do
    echo "Processing dataset: $DS..."
    
    # Define the directory and input file for this specific dataset
    DIR="dataset/positive/$DS"
    INPUT_BED="$DIR/positive_data_filtered.bed"

    # --- TEST 1: Right Padded ---
    awk 'BEGIN {OFS="\t"} {
        size = $3 - $2
        if (size < 1000) {
            new_start = $3 - 1000
            new_end = $3
            if (new_start < 0) new_start = 0
            $2 = new_start; $3 = new_end
        }
        print $0
    }' "$INPUT_BED" | bedtools intersect -a - -b "$CPG_FILE" -c > "$DIR/positive_data_filtered_right_padded.bed"

    # --- TEST 2: Center Padded ---
    awk 'BEGIN {OFS="\t"} {
        size = $3 - $2
        if (size < 1000) {
            center = int(($2 + $3) / 2)
            new_start = center - 500
            new_end = center + 500
            if (new_start < 0) new_start = 0
            $2 = new_start; $3 = new_end
        }
        print $0
    }' "$INPUT_BED" | bedtools intersect -a - -b "$CPG_FILE" -c > "$DIR/positive_data_filtered_center_padded.bed"

    # --- TEST 3: Left Padded ---
    awk 'BEGIN {OFS="\t"} {
        size = $3 - $2
        if (size < 1000) {
            new_start = $2
            new_end = $2 + 1000
            if (new_start < 0) new_start = 0
            $2 = new_start; $3 = new_end
        }
        print $0
    }' "$INPUT_BED" | bedtools intersect -a - -b "$CPG_FILE" -c > "$DIR/positive_data_filtered_left_padded.bed"

    # --- SELECTION: Compare all 3 and keep the winner ---
    # We pass the file paths into awk using the -v flag so it knows exactly where to look
    awk -v col="$TARGET_COL" \
        -v right_file="$DIR/positive_data_filtered_right_padded.bed" \
        -v center_file="$DIR/positive_data_filtered_center_padded.bed" \
        -F'\t' '{
        
        # Read Left Padded (fed via standard input)
        line1 = $0; val1 = $col
        
        # Read Right Padded
        getline line2 < right_file; split(line2, arr2, "\t"); val2 = arr2[col]
        
        # Read Center Padded
        getline line3 < center_file; split(line3, arr3, "\t"); val3 = arr3[col]
        
        # Find the max
        max_val = val1 + 0; best_line = line1
        if ((val2 + 0) > max_val) { max_val = val2 + 0; best_line = line2 }
        if ((val3 + 0) > max_val) { max_val = val3 + 0; best_line = line3 }
        
        print best_line
    }' "$DIR/positive_data_filtered_left_padded.bed" | \
    awk -F'\t' -v OFS='\t' '{NF--; print}' > "$DIR/positive_data_filtered_edited.bed"

    echo " -> Saved optimally padded regions to $DIR/positive_data_filtered_edited.bed"

done

echo "All datasets processed successfully!"
