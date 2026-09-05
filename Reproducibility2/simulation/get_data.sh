#!/bin/bash


# class="4"

for class in {1..4}; do
OUT_CASE="input/case$class"
OUT_CTR="input/ctr$class"

mkdir -p "$OUT_CASE"
mkdir -p "$OUT_CTR"

echo "Starting conversion for class.$class simulated files..."

for file in input/BL*WGBS*_background.1_DMR.class.$class.bsmooth; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        outfile="$OUT_CTR/${filename%.bsmooth}.cpg_report.txt"
        echo "Processing Control: $filename -> $outfile"
        awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3, $5, $6-$5, $4, "CGC"}' "$file" > "$outfile"
    fi
done

for file in input/FL*WGBS_background.1_DMR.class.$class.bsmooth; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        outfile="$OUT_CASE/${filename%.bsmooth}.cpg_report.txt"
        echo "Processing Case:    $filename -> $outfile"
        awk 'BEGIN{OFS="\t"} {print "chr"$1, $2, $3, $5, $6-$5, $4, "CGC"}' "$file" > "$outfile"
    fi
done

echo "Done! All class $class files have been converted and separated."
done
