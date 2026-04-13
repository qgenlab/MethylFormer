#!/bin/bash

for bed_file in dataset/DMRs/*/*methylKit*.bed; do
	echo $bed_file
	sort -k1,1 -k2,2n $bed_file | bedtools merge -i - > "${bed_file}.tmp"
	mv "${bed_file}.tmp" $bed_file
done

for bed_file in dataset/DMRs/*/*methylSig*.bed; do
	echo $bed_file
	sort -k1,1 -k2,2n $bed_file | bedtools merge -i - > "${bed_file}.tmp"
	mv "${bed_file}.tmp" $bed_file
done


mkdir -p results
######################### Overlap for NK B 

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


# positive overlap

mkdir -p results/B_NK/

for f in dataset/DMRs/B_NK/*bed; do
  echo "$f";

  base_f=$(basename "$f" .bed)

  bedtools intersect -a "$f" -b dataset/positive/NK_B/positive_data_filtered_edited.bed -F 0.05 -wo > "results/B_NK/$base_f.005.bed"
  echo "**********************************************************************"
  
  bedtools intersect -a "$f" -b dataset/positive/NK_B/positive_data_filtered_edited.bed -F 0.05 -wo

  echo "**********************************************************************"

  # bedtools intersect -a $f -b dataset/positive/NK_B/positive_data_filtered.bed -F 0.05 -wo > results/B_NK/$base_f.005.bed # The original script

  # bedtools intersect -a dataset/positive/NK_B/positive_data_filtered.bed -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c | awk -F'\t' -v OFS='\t' '$NF >= 3 {NF--; print}' | bedtools intersect -a $f -b - -F 0.05 -wo > results/B_NK/$base_f.005.bed # original script but minimum number of CpG

  # bedtools slop -i dataset/positive/NK_B/positive_data_filtered.bed -g $SCRIPT_DIR/hg38.chrom.sizes -b 0.20 -pct | bedtools intersect -a $f -b - -F 0.05 -wo > results/B_NK/$base_f.005.bed
  # bedtools slop -i dataset/positive/NK_B/positive_data_filtered.bed -g $SCRIPT_DIR/hg38.chrom.sizes -b 0.20 -pct | bedtools intersect -a - -b $f -f 0.05 -u > results/B_NK/$base_f.005.bed
done


# negative overlap (2 cells)

for f in dataset/DMRs/B_NK/*bed; do
  echo "$f";
  base_f=$(basename "$f" .bed)
  bedtools intersect -a $f -b dataset/positive/NK_B/both_datasets_data_filtered.bed -F 0.05 -wo > results/B_NK/$base_f.both_datasets_data_filtered.005.bed
done



# negative overlap

for f in dataset/DMRs/B_NK/*bed; do
  echo "$f";
  base_f=$(basename "$f" .bed)
  bedtools intersect -a $f -b dataset/negative/ctr_data_filtered_non_blood.bed -F 0.05 -wo > results/B_NK/$base_f.neg_ctr.non_blood.005.bed
done



######################### Overlap for NK Monocyte

mkdir -p results/NK_Monocytes/

# positive overlap

for f in dataset/DMRs/NK_Monocytes_res/*bed; do
  echo "$f";
  base_f=$(basename "$f" .bed)

  bedtools intersect -a "$f" -b dataset/positive/Monocyte_NK/positive_data_filtered_edited.bed -F 0.05 -wo > "results/NK_Monocytes/$base_f.005.bed"

  # bedtools intersect -a $f -b dataset/positive/Monocyte_NK/positive_data_filtered.bed -F 0.05 -wo > results/NK_Monocytes/$base_f.005.bed # The original script

  # bedtools intersect -a dataset/positive/Monocyte_NK/positive_data_filtered.bed -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c | awk -F'\t' -v OFS='\t' '$NF >= 3 {NF--; print}' | bedtools intersect -a $f -b - -F 0.05 -wo > results/NK_Monocytes/$base_f.005.bed  # original script but minimum number of CpG

  # bedtools slop -i dataset/positive/Monocyte_NK/positive_data_filtered.bed -g $SCRIPT_DIR/hg38.chrom.sizes -b 0.20 -pct | bedtools intersect -a $f -b - -F 0.05 -wo > results/NK_Monocytes/$base_f.005.bed
  # bedtools slop -i dataset/positive/Monocyte_NK/positive_data_filtered.bed -g $SCRIPT_DIR/hg38.chrom.sizes -b 0.20 -pct | bedtools intersect -a - -b $f -f 0.05 -u > results/NK_Monocytes/$base_f.005.bed
done


# negative overlap (2 cells)

for f in dataset/DMRs/NK_Monocytes_res/*bed; do
  echo "$f";
  base_f=$(basename "$f" .bed)
  bedtools intersect -a $f -b dataset/positive/Monocyte_NK/both_datasets_data_filtered.bed -F 0.05 -wo > results/NK_Monocytes/$base_f.both_datasets_data_filtered.005.bed
done



# negative overlap

for f in dataset/DMRs/NK_Monocytes_res/*bed; do
  echo "$f";
  base_f=$(basename "$f" .bed)
  bedtools intersect -a $f -b dataset/negative/ctr_data_filtered_non_blood.bed -F 0.05 -wo > results/NK_Monocytes/$base_f.neg_ctr.non_blood.005.bed
done


######################### Overlap for B Monocyte

mkdir -p results/B_Monocytes/

# positive overlap

for f in dataset/DMRs/B_Monocytes_res/*bed; do
  echo "$f";

  base_f=$(basename "$f" .bed)
 
  bedtools intersect -a "$f" -b dataset/positive/Monocyte_B/positive_data_filtered_edited.bed -F 0.05 -wo > "results/B_Monocytes/$base_f.005.bed"

  # bedtools intersect -a $f -b dataset/positive/Monocyte_B/positive_data_filtered.bed -F 0.05 -wo > results/B_Monocytes/$base_f.005.bed # The original script
  
 # bedtools intersect -a dataset/positive/Monocyte_B/positive_data_filtered.bed -b /mnt/analysis/derbelh/CpG_Without_Strand.bed -c | awk -F'\t' -v OFS='\t' '$NF >= 3 {NF--; print}' | bedtools intersect -a $f -b - -F 0.05 -wo > results/B_Monocytes/$base_f.005.bed  # original script but minimum number of CpG

  # bedtools slop -i dataset/positive/Monocyte_B/positive_data_filtered.bed -g $SCRIPT_DIR/hg38.chrom.sizes -b 0.20 -pct | bedtools intersect -a $f -b - -F 0.05 -wo > results/B_Monocytes/$base_f.005.bed
  # bedtools slop -i dataset/positive/Monocyte_B/positive_data_filtered.bed -g $SCRIPT_DIR/hg38.chrom.sizes -b 0.20 -pct | bedtools intersect -a - -b $f -f 0.05 -u > results/B_Monocytes/$base_f.005.bed
done


# negative overlap (2 cells)

for f in dataset/DMRs/B_Monocytes_res/*bed; do
  echo "$f";
  base_f=$(basename "$f" .bed)
  bedtools intersect -a $f -b dataset/positive/Monocyte_B/both_datasets_data_filtered.bed -F 0.05 -wo > results/B_Monocytes/$base_f.both_datasets_data_filtered.005.bed
done



# negative overlap

for f in dataset/DMRs/B_Monocytes_res/*bed; do
  echo "$f";
  base_f=$(basename "$f" .bed)
  bedtools intersect -a $f -b dataset/negative/ctr_data_filtered_non_blood.bed -F 0.05 -wo > results/B_Monocytes/$base_f.neg_ctr.non_blood.005.bed
done

