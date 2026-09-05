

project_path=$1

paths=("./B_Monocytes_res/" "./NK_Monocytes_res/" "./B_NK/" "./KOB_WTA/" "./AD_res/")



for dataset in "${paths[@]}"; do
        echo "Processing " $dataset
        ./prepare_benchmark.sh $project_path/$dataset/methylkit/new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.bed \
                               $project_path/$dataset/methylsig/new2.methylSig..hg38.window.1000.bed \
                               $project_path/$dataset/dss/new2.dss.CpG.hg38DMR.bed \
                               $project_path/$dataset/bsseq/new2.dsseq..hg38.DMR.bed \
                               $project_path/$dataset/DiffMethylTools/data/generate_DMR_0.bed \
                               $project_path/$dataset/DiffMethylTools_dl/data/filtered_bed/Gastric_cancer_l1_pen_0_1_scaled_sum_1_2.bed \
                               $project_path/$dataset/benchmark_3
done

