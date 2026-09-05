#!/bin/bash

if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is required but not installed or not in PATH." >&2
    exit 1
fi

paths=(
    "./B_Monocytes_res"
    "./NK_Monocytes_res"
    "./B_NK"
    "./KOB_WTA"
    "./AD_res"
)

benchmarks=(
    "benchmark_3.new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.benchmark.bed"
    "benchmark_3.new2.methylSig..hg38.window.1000.benchmark.bed"
    "benchmark_3.new2.dss.CpG.hg38DMR.benchmark.bed"
    "benchmark_3.new2.dsseq..hg38.DMR.benchmark.bed"
    "benchmark_3.generate_DMR_0.benchmark.bed"
    "benchmark_3.generate_DMR_CPD.benchmark.bed"
)

results=(
    "methylkit/new2.methylkit..destranded.CpG.hg38.window.1000.step.500.cov.5.bed"
    "methylsig/new2.methylSig..hg38.window.1000.bed"
    "dss/new2.dss.CpG.hg38DMR.bed"
    "bsseq/new2.dsseq..hg38.DMR.bed"
    "DiffMethylTools/data/generate_DMR_0.bed"
    "DiffMethylTools_dl/data/generate_DMR_CPD.bed"
)

echo -e "Path\tBenchmark_File\tResult_File\tTotal_Benchmark_Regions\tTotal_Result_Regions\tOverlapping_Result_Regions\tNon_Overlapping_Result_Regions"

for p in "${paths[@]}"; do
    echo $p
    clean_p="${p%/}"
    for i in "${!benchmarks[@]}"; do
        bench_file="${clean_p}/${benchmarks[$i]}"
        res_file="${clean_p}/${results[$i]}"
        if [[ -f "$bench_file" && -f "$res_file" ]]; then
            total_bench=$(wc -l < "$bench_file" | tr -d ' ')
            total_res=$(wc -l < "$res_file" | tr -d ' ')
            overlap_res=$(bedtools intersect -a "$res_file" -b "$bench_file" -f 0.1 -u | wc -l | tr -d ' ')
            non_overlap_res=$(bedtools intersect -a "$res_file" -b "$bench_file" -v | wc -l | tr -d ' ')
            echo -e "${clean_p}\t${benchmarks[$i]}\t${results[$i]}\t${total_bench}\t${total_res}\t${overlap_res}\t${non_overlap_res}"
        else
            echo -e "${clean_p}\t${benchmarks[$i]}\t${results[$i]}\tFILE_NOT_FOUND\tFILE_NOT_FOUND\tNA\tNA"
        fi
    done
done

