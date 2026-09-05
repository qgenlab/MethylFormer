#!/bin/bash

set -e

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 file1.bed file2.bed file3.bed file4.bed file5.bed file6.bed output"
    exit 1
fi

TMP_DIR=$(mktemp -d)
echo "Created temporary directory for preprocessing: $TMP_DIR"

PROCESSED_FILES=()

for i in {1..6}; do
	INFILE="${!i}"
	filename="${INFILE##*/}"
	clean_name="${filename%.*}"
	OUTFILE="$TMP_DIR/${clean_name}.bed"
	echo "Preprocessing file $i: $INFILE -> $OUTFILE ..."
	bedtools sort -i "$INFILE" | bedtools merge -i - -d 10 > "$OUTFILE"
	PROCESSED_FILES+=("$OUTFILE")
done


# for i in {1..6}; do
#     INFILE="${!i}"
#     OUTFILE="$TMP_DIR/processed_file_${i}.bed"
#    echo "Preprocessing file $i: $INFILE ..."
#    bedtools sort -i "$INFILE" | \
#        bedtools merge -i - -d 10 > "$OUTFILE"
#    PROCESSED_FILES+=("$OUTFILE")
# done

echo "--------------------------------------------------"
echo "Finding regions present in at least 4 files..."


total_items=${#PROCESSED_FILES[@]}

for ((i=0; i<total_items; i++)); do
	current_file="${PROCESSED_FILES[0]}"
	PROCESSED_FILES=("${PROCESSED_FILES[@]:1}")
	echo "Processing: $current_file"
	echo "List to be processed: ${PROCESSED_FILES[*]}"
	filename="${current_file##*/}"
	basename="${filename%.*}"
	bedtools multiinter -i "${PROCESSED_FILES[@]}" | \
    		awk '$4 >= 3 {print $1"\t"$2"\t"$3}' | bedtools merge -i - -d 100 > $7.$basename.benchmark.bed
	PROCESSED_FILES+=("$current_file")
done


rm -rf "$TMP_DIR"

echo "Done! Final merged output saved to: $7"
