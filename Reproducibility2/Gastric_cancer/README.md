## Download data
To download data run get_data.sh.

## bismark_pipe.sh Command

bismark_pipe.sh input_folder output_folder threads fastqc_env java_env bismark_env sambamba_env ref_path
**Example**:
To run the command in the current directory with 30 threads:
`../bismark_pipe.sh . . 30 lrst_py39-1 java lrst_py39-1 sambamba /home/bismark/hg38/`
Here, lrst_py39-1, java, lrst_py39-1, and sambamba represent the required environments.

**Processing Steps**:
Step 0 (One-time setup):
Run bismark_genome_preparation on the reference genome using Bowtie2.

Note: "no_rand" indicates that random chromosomes have been filtered out from the reference genome.
Step 1:
Merge all input files within the specified folder into a single .fq file.
Note: The script assumes paired-end sequencing data, where paired reads end with _1 and _2.

Step 2:
Perform quality control checks on the input data using FastQC.

Step 3:
Trim the reads using Trimmomatic with Phred+33 quality encoding. The trimming criteria are:
- Remove bases at the end of reads with a quality score below 20.
- Discard reads shorter than 30 bases after trimming.
Retain both paired reads (saved as _paired.fq) and unpaired reads (saved as _unpaired.fq).

Step 4:
Run FastQC again on the trimmed data to assess quality.

Step 5:
Align the reads using Bismark with Bowtie2 as the alignment algorithm.

Step 6:
Identify and mark duplicate reads in the BAM file generated from Step 5.

Step 7:
Filter aligned reads in the BAM file using SAMtools, retaining only reads with a mapping quality ≥ 10.
The -F flag is used to exclude reads with the following flags:

```
4    -> Unmapped reads
8    -> Mate is unmapped
256  -> Secondary alignments
512  -> Fails quality control
1024 -> PCR/Optical duplicates
```

Step 8:
Extract methylation information from the Bismark-aligned BAM file using bismark_methylation_extractor.

Generates a detailed methylation report including all cytosine contexts (CpG, CHG, CHH).

