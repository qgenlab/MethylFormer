

if [ "$#" -ne 8 ]; then
	echo "Step 0 (One-time setup):"
	echo "Run bismark_genome_preparation on the reference genome using Bowtie2. The output folder will be your 8th parameter."
	echo "Step 1"
	echo "To run the command in the current directory with 30 threads:"
	echo "bismark_pipe.sh input_folder_path output_folder_path nbr_of_threads fastqc_env java_env bismark_env sambamba_env"
	exit 0
fi

input=$1
output=$2
threads=$3
env1=$4
env2=$5
env3=$6
env4=$7
ref=$8

# ./get_data.sh $inf $sup $output > $output/fastq.log

cat `ls $input/*_1.fastq` > $output/reads_1.fq&
cat `ls $input/*_2.fastq` > $output/reads_2.fq&

wait

eval "$(conda shell.bash hook)"
conda activate $env1


fastqc -t 10 $output/reads_1.fq&
fastqc -t 10 $output/reads_2.fq&


wc $output/*.fq >> $output/mreadme.wc.txt&
wc $output/*.fastq >> $output/mreadme.wc.txt&

eval "$(conda shell.bash hook)"
conda activate $env2

mkdir -p $output/trimmed

java -jar ~/bin/trimmomatic.jar PE -threads $threads -phred33 $output/reads_1.fq $output/reads_2.fq $output/trimmed/output_forward_paired.fq $output/trimmed/output_forward_unpaired.fq $output/trimmed/output_reverse_paired.fq $output/trimmed/output_reverse_unpaired.fq TRAILING:20 MINLEN:30 > trimmomatic.log.txt 2> trimmomatic.err.log


eval "$(conda shell.bash hook)"
conda activate $env1

fastqc -t 10 $output/trimmed/output_forward_paired.fq&
fastqc -t 10 $output/trimmed/output_reverse_paired.fq&

eval "$(conda shell.bash hook)"
conda activate $env3

mkdir -p $output/trimmed/tmp

/mnt/analysis/derbelh/bismark/Bismark-0.24.2/bismark --bowtie2 $ref -1 $output/trimmed/output_forward_paired.fq -2 $output/trimmed/output_reverse_paired.fq --parallel $threads -o $output/trimmed/ --temp_dir $output/trimmed/tmp/ >> $output/trimmed/bismark.log 2> $output/trimmed/bismark.err.log

eval "$(conda shell.bash hook)"
conda activate $env4

sambamba markdup -l 1 -t $threads --sort-buffer-size 16000 --overflow-list-size 10000000 $output/trimmed/output_forward_paired_bismark_bt2_pe.bam $output/trimmed/marked_output_forward_paired_bismark_bt2_pe.bam

eval "$(conda shell.bash hook)"
conda activate $env3


samtools view -q 10 -F 1796 $output/trimmed/marked_output_forward_paired_bismark_bt2_pe.bam -o $output/trimmed/filtered_marked_output_forward_paired_bismark_bt2_pe.bam

/mnt/analysis/derbelh/bismark/Bismark-0.24.2/bismark_methylation_extractor -p --comprehensive --no_overlap --cytosine_report --genome_folder $ref --zero_based $output/trimmed/filtered_marked_output_forward_paired_bismark_bt2_pe.bam --parallel $threads -o $output/trimmed/>> $output/trimmed/bismark_methylation_extractor.log 2> $output/trimmed/bismark_methylation_extractor.err.log
