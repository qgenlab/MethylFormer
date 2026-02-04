

# example usage: ./softlink.sh "/Blood-NK/*/data_CpG_report.txt" test_nk_mono/case
FILES=$1
OUTPUT=$2

for file in $FILES; do
    [ -e "$file" ] || continue
    sample_id=$(basename "$(dirname "$file")")
    filename=$(basename "$file")
    new_name="${sample_id}_${filename}"
    echo $OUTPUT
    echo "$OUTPUT/$new_name"
    ln -sf "$file" "$OUTPUT/$new_name"
    
    echo "Linked: $new_name"
done
