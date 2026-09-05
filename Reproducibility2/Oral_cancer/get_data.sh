

wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542054/suppl/GSM6542054%5FOC2%5FMergedCG%5FmC%5F10x%2Etxt%2Egz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542050/suppl/GSM6542050%5FOC1%5FMergedCG%5FmC%5F10x%2Etxt%2Egz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542058/suppl/GSM6542058%5FOC3%5FMergedCG%5FmC%5F10x%2Etxt%2Egz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542062/suppl/GSM6542062%5FOC4%5FMergedCG%5FmC%5F10x%2Etxt%2Egz -P ./case/



wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542059/suppl/GSM6542059%5FON3%5FMergedCG%5FmC%5F10x%2Etxt%2Egz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542063/suppl/GSM6542063%5FON4%5FMergedCG%5FmC%5F10x%2Etxt%2Egz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542055/suppl/GSM6542055%5FON2%5FMergedCG%5FmC%5F10x%2Etxt%2Egz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6542nnn/GSM6542051/suppl/GSM6542051%5FON1%5FMergedCG%5FmC%5F10x%2Etxt%2Egz -P ./ctr/



gunzip ./case/*
gunzip ./ctr/*


for f in */*.txt; do
    echo $f
    ./preprocess.sh $f
done


echo "Downloading Benchmark..."

wget https://methmarkerdb.hzau.edu.cn/memarkerdata/result_0808/Oral_cancer/Oral_cancer.SMART.dmr.merge.DMRanno.txt


awk 'NR > 1 && $1 !~ /_random/ {print $1 "\t" $2 "\t" $3 "\t" $12}' Oral_cancer.SMART.dmr.merge.DMRanno.txt | bedtools sort -i - > Oral_cancer.SMART.dmr.merge.DMRanno.bed
