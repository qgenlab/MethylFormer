

# path="/mnt/analysis/derbelh/github_test/test_used_for_paper_final/simulation4/"

path=$1


dmr_files=(
    "$path/dss/new2.dss.CpG.hg38DMR.csv"
    "$path/methylkit/new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.csv"
    "$path/methylsig/new2.methylSig..hg38.window.1000.csv"
    "$path/bsseq/new2.dsseq..hg38.DMR.csv"
    "$path/DiffMethylTools/data/generate_DMR_0.csv" # old dl clustring algorithm
    "$path/DiffMethylTools_dl/data/generate_DMR_CPD_0.csv"
)

preprocessed_dmr_files=(
)

#dmr_file="/home/derbelh/analysis/data_test/metilene_simulation/metilene_DMR_simulation_scripts/simulationdata/DMRs/DMR.artif.chr10.WGBS_background.1_DMR.class.4.bed"

dmr_file=$2


for file in "${dmr_files[@]}"; do
    PATH_FILE=$(dirname $file)
    FILE_NAME=$(basename "$file")
    NAME_ONLY="${FILE_NAME%.*}"
    echo $NAME_ONLY
    # awk -v OFS="\t" -F, 'NR > 1 {print $1, $2, $3}' $file | grep -v -E 'random|Un|alt' | sed 's/"//g' | bedtools sort -i - > $PATH_FILE/$NAME_ONLY.preprocessed.bed
    preprocessed_dmr_files+=("$PATH_FILE/$NAME_ONLY.bed")
done



echo "testing with -f 0.1"

for file in "${preprocessed_dmr_files[@]}"; do
    echo "Processing: $file"
    echo "Numerator for recall"
    bedtools intersect -a $dmr_file -b $file -f 0.1 -u | wc -l

    echo "denominator for recall"
    wc -l "$dmr_file"

    echo "Numerator for precision"
    bedtools intersect -a $file -b $dmr_file -f 0.1 -u | wc -l
    
    echo "denominator for precision"
    wc -l $file

done


echo "testing with -f 0.25"

for file in "${preprocessed_dmr_files[@]}"; do
    echo "Processing: $file"
    echo "Numerator for recall"
    bedtools intersect -a $dmr_file -b $file -f 0.25 -u | wc -l

    echo "denominator for recall"
    wc -l "$dmr_file"

    echo "Numerator for precision"
    bedtools intersect -a $file -b $dmr_file -f 0.25 -u | wc -l

    echo "denominator for precision"
    wc -l $file

done




echo "testing with -f 0.5"

for file in "${preprocessed_dmr_files[@]}"; do
    echo "Processing: $file"
    echo "Numerator for recall"
    bedtools intersect -a $dmr_file -b $file -f 0.5 -u | wc -l

    echo "denominator for recall"
    wc -l "$dmr_file"

    echo "Numerator for precision"
    bedtools intersect -a $file -b $dmr_file -f 0.5 -u | wc -l

    echo "denominator for precision"
    wc -l $file
done


