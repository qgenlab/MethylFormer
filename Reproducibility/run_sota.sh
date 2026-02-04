#!/bin/bash


Help()
{
   # Display Help
   echo "This script generates the results for all tools."
   echo
   echo "Syntax: ./run.sh [-t|i|f|o|r|a|b|s]"
   echo "options:"
   echo "t     Number of threads. [default=4]"
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

mkdir -p $output
mkdir -p $output/case
mkdir -p $output/ctr


echo "Converting files..."
case=$case/*
ctr=$ctr/*

for f in $case; do
	python ./convert/MethylConverter.py $f $input_format $output/case;
done

for f in $ctr; do
        python ./convert/MethylConverter.py $f $input_format $output/ctr;
done


echo "running DSS..."
mkdir -p $output/dss
ctr_list=$(ls $output/ctr/*converted.dss.txt | paste -sd "," -)
case_list=$(ls $output/case/*converted.dss.txt | paste -sd "," -)

Rscript run/run_dss.r --case "$case_list" --control "$ctr_list" -d TRUE --assembly hg38 --output $output/dss/


eval "$(conda shell.bash hook)"
conda activate mekit

echo "running Methylkit..."
if [[ $input_format == "CR" ]]; then
	ctr_list=$(ls $ctr/* | paste -sd "," -)
	case_list=$(ls $case/* | paste -sd "," -)
else
        ctr_list=$(ls $output/ctr/*converted.cytosine_report.txt | paste -sd "," -)
        case_list=$(ls $output/case/*converted.cytosine_report.txt | paste -sd "," -)
fi

mkdir -p $output/methylkit/
Rscript run/run_methylkit.r --case "$case_list" --control "$ctr_list" -w FALSE --destrand TRUE --assembly hg38 --context CpG -m 2 --output $output/methylkit/ -b 1

Rscript run/run_methylkit.r --case "$case_list" --control "$ctr_list" -w TRUE --destrand TRUE --assembly hg38 --context CpG -m 2 --output $output/methylkit/ -b 1


eval "$(conda shell.bash hook)"
conda activate r441


echo "running MethylSig..."
if [[ $input_format == "CR" ]]; then
        ctr_list=$(ls $ctr/* | paste -sd "," -)
        case_list=$(ls $case/* | paste -sd "," -)
else
        ctr_list=$(ls $output/ctr/*converted.cytosine_report.txt | paste -sd "," -)
        case_list=$(ls $output/case/*converted.cytosine_report.txt | paste -sd "," -)

fi

mkdir -p $output/methylsig/
Rscript run/run_methylSig.r --case "$case_list" --control "$ctr_list" --cov 10 --spcov1 1 --spcov2 1 -w FALSE --output $output/methylsig/
Rscript run/run_methylSig.r --case "$case_list" --control "$ctr_list" --cov 10 --spcov1 1 --spcov2 1 -w TRUE --output $output/methylsig/

echo "running bsseq..."
if [[ $input_format == "CR" ]]; then
        ctr_list=$(ls $ctr/* | paste -sd "," -)
        case_list=$(ls $case/* | paste -sd "," -)
else
        ctr_list=$(ls $output/ctr/*converted.cytosine_report.txt | paste -sd "," -)
        case_list=$(ls $output/case/*converted.cytosine_report.txt | paste -sd "," -)

fi

mkdir -p $output/bsseq/
Rscript run/run_bsseq.r --case "$case_list" --control "$ctr_list" --spcov1 2 --spcov2 2  --cov1 10 --cov2 10 -w FALSE --output $output/bsseq/



echo "running DiffMethylTools..."

if [[ $input_format == "BED" ]]; then
        ctr_list=$(ls $ctr/* | paste -sd " " -)
        case_list=$(ls $case/* | paste -sd " " -)
else
        ctr_list=$(ls $output/ctr/*_converted.bed | paste -sd " " -)
        case_list=$(ls $output/case/*_converted.bed | paste -sd " " -)

fi

mkdir -p $output/DiffMethylTools/
python ../DiffMethylTools.py all_analysis --case "$case_list" --control "$ctr_list" --input_format BED --ref_folder hg38




###########################################################
###########################################################
###########################################################
#
#if [[ $input_format == "BED" ]]; then
#	case=$case/*.bed
#	ctr=$ctr/*.bed
#	echo "Converting from BED to DSS format..."
#	for f in $case; do
#		./convert/bed2dss.py $f $output/case
#	done
#        for f in $ctr; do
#                ./convert/bed2dss.py $f $output/ctr
#        done
#
#fi
#	
#	
#	
#if [[ $input_format == "CR" ]]; then	
#	echo "Converting from  format..."
#        for f in $case; do
#                ./convert/bed2dss.py $f $output/case
#        done
#        for f in $ctr; do
#                ./convert/bed2dss.py $f $output/ctr
#        done
#fi
#
#
#if [[ $input_format == "MB" ]]; then
#	echo "Converting from format..."
#        for f in $case; do
#                ./convert/bed2dss.py $f $output/case
#        done
#        for f in $ctr; do
#                ./convert/bed2dss.py $f $output/ctr
#        done
#fi
#
