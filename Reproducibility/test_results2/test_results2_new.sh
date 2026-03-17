
# /mnt/analysis/derbelh/github_test/DiffMethylTools2/Reproducibility/test_results2/test_results2.sh /mnt/analysis/derbelh/github_test/test1/B_NK/ /mnt/analysis/derbelh/github_test/test1/NK_Monocytes_res/ /mnt/analysis/derbelh/github_test/test1/B_Monocytes_res/

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

echo "Running benchmark generation..."
python "$SCRIPT_DIR/generate_benchmark.py" "$1" "$2" "$3"


eval "$(conda shell.bash hook)"
conda activate lrst_py39-1

echo "Running edit benchmark..."
"$SCRIPT_DIR/edit_benchmark2.sh"

echo "Running overlap..."
"$SCRIPT_DIR/overlap2.sh"

eval "$(conda shell.bash hook)"
conda activate base


echo "Running plotting..."
python "$SCRIPT_DIR/plot.py"

echo "Pipeline complete!"
