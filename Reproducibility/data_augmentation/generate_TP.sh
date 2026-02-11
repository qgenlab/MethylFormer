all_merged_2=(
    "$HOME/analysis/DNA_data_pre/pipeline_project/intersect_gt/2/dss._DiffMethylTools..bed" 
    "$HOME/analysis/DNA_data_pre/pipeline_project/intersect_gt/2/DiffMethylTools._dsseq..bed" 
    "$HOME/analysis/DNA_data_pre/pipeline_project/intersect_gt/2/dss._dsseq..bed" 
    "$HOME/analysis/DNA_data_pre/pipeline_project/intersect_gt/2/methylkit._dsseq..bed" 
    "$HOME/analysis/DNA_data_pre/pipeline_project/intersect_gt/2/methylkit._DiffMethylTools..bed" 
    "$HOME/analysis/DNA_data_pre/pipeline_project/intersect_gt/2/methylSig._DiffMethylTools..bed" 
    "$HOME/analysis/DNA_data_pre/pipeline_project/intersect_gt/2/methylSig._dsseq..bed" 
)


CpG="/mnt/analysis/derbelh/CpG_Without_Strand.bed"

CpG_output_files=()

for f in "${all_merged_2[@]}"; do
    echo $f
    filename=$(basename "$f")
    basename_no_ext="${filename%.*}"
    bedtools intersect -a $CpG -b $f > "./"$basename_no_ext"_pos.bed"
    CpG_output_files+=("./"$basename_no_ext"_pos.bed")
done

for f in "${CpG_output_files[@]}"; do
    tmp="${f}.tmp"
    bedtools sort -i "$f" > "$tmp" && mv "$tmp" "$f"
done


bedtools multiinter -i "${CpG_output_files[@]}" > union_intersect_2.bed


bedtools merge -d 100 -c 1 -o count -i union_intersect_2.bed > union_intersect_2.bed_with_counts.bed

awk '($3 - $2) >= 50 && $4 >= 3' union_intersect_2.bed_with_counts.bed > union_intersect_2.bed_with_counts_filtered.bed

awk '{$2 = ($2 - 1000 < 0 ? 0 : $2 - 1000); $3 = $3 + 1000; print}' OFS='\t' union_intersect_2.bed_with_counts_filtered.bed |  sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o sum > union_intersect_2.bed_with_counts_filtered_1kbp.bed

