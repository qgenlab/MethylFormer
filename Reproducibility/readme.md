# State-of-the-Art (SOTA) Comparison Pipeline
This repository provides an automated pipeline to run and evaluate state-of-the-art tools for differential methylation analysis plus the proposed MethylFormer model.

## 1. Data Preparation (Optional)
To keep your workspace clean and avoid copying large datasets, it is highly recommended to create symbolic links to your input data before running the pipeline.

You can use the provided `softlink.sh` script to quickly link your files into your working directory:
```
./softlink.sh "path/to/case/files/*.txt" path/to/output_folder/
```

## 2. Generate the Results
The `run_sota.sh` script is the core of the pipeline. It takes your raw data and automatically converts it into the specific input formats required by each individual tool before executing them.

Usage:
```
./run_sota.sh -i <input_format> -a <case_folder> -r <control_folder> -o <output_folder/>
Parameters:
```

-i : The format of your input data. Must be one of the following:

BED — Standard BED format

CR — Bismark CpG Report

MB — methylBase (methylKit) format

-a : Path to the folder containing the Case files.

-r : Path to the folder containing the Control files.

-o : Path to the desired output directory.

## 3. Analyze the Performance
Once the SOTA tools have finished running, use the `analyse_results.py` script to evaluate their performance.

This script will automatically detect the .yml configuration file generated in the previous step, locate all the predicted results, and calculate the final performance metrics.

Usage:
```
python analyse_results.py path/to/output_folder/
```
Note: The path provided here should be the exact same output_folder/ you specified in Step 2. The script will output comprehensive performance tables based on the tools' predictions.
