INPUT_FILES="/mnt/analysis/derbelh/DNA_data_pre/MO3_ND3_13.02.24/*/filtered/filtered_CpG/new_destranded_*_5mc.new.methyl1_filtered_1.CpG.bed"


mkdir -p tmp

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -P tmp/
gunzip tmp/hg38.fa.gz


python generate_wg_cpg_sorted_wostrand.py tmp/hg38.fa tmp/CpG_Without_Strand.bed
python generate_regions.py tmp/CpG_Without_Strand.bed

#
##
all_merged_2=(
    "dmrs/dss._DiffMethylTools..bed" 
    "dmrs/DiffMethylTools._dsseq..bed" 
    "dmrs/dss._dsseq..bed" 
    "dmrs/methylkit._dsseq..bed" 
    "dmrs/methylkit._DiffMethylTools..bed" 
    "dmrs/methylSig._DiffMethylTools..bed" 
    "dmrs/methylSig._dsseq..bed" 
)
##
##
CpG="tmp/CpG_Without_Strand.bed"
##
CpG_output_files=()
##
for f in "${all_merged_2[@]}"; do
    echo $f
    filename=$(basename "$f")
    basename_no_ext="${filename%.*}"
    bedtools intersect -a $CpG -b $f > "tmp/"$basename_no_ext"_pos.bed"
    CpG_output_files+=("tmp/"$basename_no_ext"_pos.bed")
done
##
for f in "${CpG_output_files[@]}"; do
    tmp="${f}.tmp"
    bedtools sort -i "$f" > "$tmp" && mv "$tmp" "$f"
done
##
##
bedtools multiinter -i "${CpG_output_files[@]}" > tmp/union_intersect_2.bed
##
##
bedtools merge -d 100 -c 1 -o count -i tmp/union_intersect_2.bed > tmp/union_intersect_2.bed_with_counts.bed
##
awk '($3 - $2) >= 50 && $4 >= 3' tmp/union_intersect_2.bed_with_counts.bed > tmp/union_intersect_2.bed_with_counts_filtered.bed
##
awk '{$2 = ($2 - 1000 < 0 ? 0 : $2 - 1000); $3 = $3 + 1000; print}' OFS='\t' tmp/union_intersect_2.bed_with_counts_filtered.bed |  sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o sum > tmp/union_intersect_2.bed_with_counts_filtered_1kbp.bed
##
python split_regions.py
#

for f in $INPUT_FILES; do
	echo $f
	python simulation_1.py $f
done


python ./map_reg_to_pos.py

mkdir -p ctr
mkdir -p rev_ctr


for f in $INPUT_FILES; do
	echo $f
	python generate_ctr.py $f
done


mkdir -p case


python ./smoothed.py


echo "start generating outputs"
python generate_output.py

echo "start merging df"

python merge_out_meth.py

echo "start smoothing outputs"
python smooth_cohn_d.py
# python map_reg_to_pos.py #output of simulation_1.py

