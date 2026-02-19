#!/bin/bash


Help()
{
   # Display Help
   echo "This script generates the results for all tools."
   echo
   echo "Syntax: ./run.sh [-t|i|f|o|r|a|b|s]"
   echo "options:"
   echo "t     Number of threads. [Not yet implemented]"
   echo "i     The input file format. [PS: BED for bed format, CR for Bismark CpG Report or MB for methylBase (methylKit) format]"
   echo "a     The input case folder."
   echo "r     The input ctr folder."
   echo "o     The output folder."
   echo
}


while getopts "ht:i:a:r:o:" opt; do
		    case $opt in
		    h) Help
		    exit;;
		    t) threads=${OPTARG};;
		    i) input_format=${OPTARG};;
		    a) case=${OPTARG};;
		    r) ctr=${OPTARG};;
		    o) output=${OPTARG};;
		    esac
done

echo $threads $input_format $case $ctr $output



SCRIPT_DIR=$(dirname "$0")

mkdir -p "$output"
mkdir -p "$output/case"
mkdir -p "$output/ctr"


echo "" > $output/config.yml


#
#
echo "Converting files..."
case=$case/*
ctr=$ctr/*
#
for f in $case; do
	python $SCRIPT_DIR/convert/MethylConverter.py $f $input_format $output/case;
done

for f in $ctr; do
        python $SCRIPT_DIR/convert/MethylConverter.py $f $input_format $output/ctr;
done

eval "$(conda shell.bash hook)"
conda activate r441


echo "running DSS..."
mkdir -p $output/dss
ctr_list=$(ls $output/ctr/*converted.dss.txt | paste -sd "," -)
case_list=$(ls $output/case/*converted.dss.txt | paste -sd "," -)

Rscript $SCRIPT_DIR/run/run_dss.r --case "$case_list" --control "$ctr_list" -d TRUE --assembly hg38 --output $output/dss/


echo "DSS_dmr: "$output/dss/new2.dss.CpG.hg38DMR.csv >> $output/config.yml
echo "DSS_dml: "$output/dss/new2.dss.CpG.hg38DML.csv >> $output/config.yml


eval "$(conda shell.bash hook)"
conda activate mekit

echo "running Methylkit..."
if [[ $input_format == "CR" ]]; then
	ctr_list=$(ls $ctr | paste -sd "," -)
	case_list=$(ls $case | paste -sd "," -)
else
        ctr_list=$(ls $output/ctr/*converted.cytosine_report.txt | paste -sd "," -)
        case_list=$(ls $output/case/*converted.cytosine_report.txt | paste -sd "," -)
fi

mkdir -p $output/methylkit/
Rscript $SCRIPT_DIR/run/run_methylkit.r --case "$case_list" --control "$ctr_list" -w FALSE --destrand TRUE --assembly hg38 --context CpG -m 2 --output $output/methylkit/ -b 1

Rscript $SCRIPT_DIR/run/run_methylkit.r --case "$case_list" --control "$ctr_list" -w TRUE --destrand TRUE --assembly hg38 --context CpG -m 2 --output $output/methylkit/ -b 1

echo "methylKit_dmr: "$output/methylkit/new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.10.csv >> $output/config.yml
echo "methylKit_dml: "$output/methylkit/new2.methylkit..destranded.CpG.hg38.non.window.cov.10.csv  >> $output/config.yml


eval "$(conda shell.bash hook)"
conda activate r441


echo "running MethylSig..."
if [[ $input_format == "CR" ]]; then
        ctr_list=$(ls $ctr | paste -sd "," -)
        case_list=$(ls $case | paste -sd "," -)
else
        ctr_list=$(ls $output/ctr/*converted.cytosine_report.txt | paste -sd "," -)
        case_list=$(ls $output/case/*converted.cytosine_report.txt | paste -sd "," -)

fi

mkdir -p $output/methylsig/
Rscript $SCRIPT_DIR/run/run_methylSig.r --case "$case_list" --control "$ctr_list" --cov 10 --spcov1 1 --spcov2 1 -w FALSE --output $output/methylsig/
Rscript $SCRIPT_DIR/run/run_methylSig.r --case "$case_list" --control "$ctr_list" --cov 10 --spcov1 1 --spcov2 1 -w TRUE --output $output/methylsig/


echo "methylSig_dml: "$output/methylsig/new2.methylSig..hg38.non.window.csv >> $output/config.yml
echo "methylSig_dmr: "$output/methylsig/new2.methylSig..hg38.window.1000.csv >> $output/config.yml



echo "running bsseq..."
if [[ $input_format == "CR" ]]; then
        ctr_list=$(ls $ctr | paste -sd "," -)
        case_list=$(ls $case | paste -sd "," -)
else
        ctr_list=$(ls $output/ctr/*converted.cytosine_report.txt | paste -sd "," -)
        case_list=$(ls $output/case/*converted.cytosine_report.txt | paste -sd "," -)

fi

mkdir -p $output/bsseq/
Rscript $SCRIPT_DIR/run/run_bsseq.r --case "$case_list" --control "$ctr_list" --spcov1 2 --spcov2 2  --cov1 10 --cov2 10 --output $output/bsseq/

echo "BSseq: "$output/bsseq/new2.dsseq..hg38.DMR.csv >> $output/config.yml


eval "$(conda shell.bash hook)"
conda activate /mnt/analysis/derbelh/.local/share/mamba/envs/DiffMethylTools
echo "running DiffMethylTools..."

if [[ $input_format == "BED" ]]; then
        ctr_list=$(ls $ctr | paste -sd " " -)
        case_list=$(ls $case | paste -sd " " -)
else
        ctr_list=$(ls $output/ctr/*_converted.bed | paste -sd " " -)
        case_list=$(ls $output/case/*_converted.bed | paste -sd " " -)

fi



mkdir -p $output/DiffMethylTools/
python $SCRIPT_DIR/../DiffMethylTools.py all_analysis --case_data_file $case_list --ctr_data_file $ctr_list --input_format BED --ref_folder hg38 --results_path $output/DiffMethylTools/

echo "DiffMethylTools_dmr: "$output/DiffMethylTools/data/generate_DMR_0.csv >> $output/config.yml
echo "DiffMethylTools_dml: "$output/DiffMethylTools/data/filters.csv >> $output/config.yml


echo "running DiffMethylTools with Deep learning..."

if [[ $input_format == "BED" ]]; then
        ctr_list=$(ls $ctr | paste -sd " " -)
        case_list=$(ls $case | paste -sd " " -)
else
        ctr_list=$(ls $output/ctr/*_converted.bed | paste -sd " " -)
        case_list=$(ls $output/case/*_converted.bed | paste -sd " " -)

fi

mkdir -p $output/DiffMethylTools_dl/
python $SCRIPT_DIR/../DiffMethylTools.py all_analysis_dl --case_data_file $case_list --ctr_data_file $ctr_list --input_format BED --ref_folder hg38 --model_path $SCRIPT_DIR/../bin/DL_model_state.pth --results_path $output/DiffMethylTools_dl/

# echo "dl_dmr_040: "$output/DiffMethylTools_dl/data/ >> $output/config.yml
echo "dl_dml: "$output/DiffMethylTools_dl/data/filters_dl.csv >> $output/config.yml

echo "dl_dmr_035: "$output/DiffMethylTools_dl/data/generate_DMR_0.csv >> $output/config.yml


echo "all_data: "$output/DiffMethylTools/data/merge_tables.csv >> $output/config.yml
