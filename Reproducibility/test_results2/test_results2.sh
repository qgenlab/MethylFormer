
# /mnt/analysis/derbelh/github_test/DiffMethylTools2/Reproducibility/test_results2/test_results2.sh /mnt/analysis/derbelh/github_test/test1/B_NK/ /mnt/analysis/derbelh/github_test/test1/NK_Monocytes_res/ /mnt/analysis/derbelh/github_test/test1/B_Monocytes_res/

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "Running benchmark generation..."
python "$SCRIPT_DIR/generate_benchmark.py" "$1" "$2" "$3"


eval "$(conda shell.bash hook)"
conda activate lrst_py39-1


echo "Running overlap..."
"$SCRIPT_DIR/overlap.sh"

eval "$(conda shell.bash hook)"
conda activate base



echo "Running plotting..."
#python "$SCRIPT_DIR/plot.py"
python "$SCRIPT_DIR/01_generate_benchmarks.py"

eval "$(conda shell.bash hook)"
conda activate lrst_py39-1

# "$SCRIPT_DIR/02_filter_TP.sh" # OG script
"$SCRIPT_DIR/12_filter_TP.sh"

eval "$(conda shell.bash hook)"
conda activate base

python "$SCRIPT_DIR/03_calculate_metrics.py"

echo "Pipeline complete!"
