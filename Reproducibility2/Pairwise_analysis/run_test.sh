

path_to_results_folder=$1
SCRIPT_DIR=$(dirname "$0")


python $SCRIPT_DIR/../../Reproducibility/analyse_results.py $path_to_results_folder/NK_Monocytes_res &


python $SCRIPT_DIR/../../Reproducibility/analyse_results.py $path_to_results_folder/B_Monocytes_res &


python $SCRIPT_DIR/../../Reproducibility/analyse_results.py $path_to_results_folder/B_NK &


python $SCRIPT_DIR/../../Reproducibility/analyse_results.py $path_to_results_folder/KOB_WTA &


python $SCRIPT_DIR/../../Reproducibility/analyse_results.py $path_to_results_folder/AD_res &


wait
