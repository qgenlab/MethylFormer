
filepath=$1

folder="${filepath%/*}"
filename="${filepath##*/}"
gsm_id="${filename%%_*}"


awk -v OFS='\t' 'NR > 1 {print $1, $2, $3, $7, $6 - $7, $4, "CGC"}' $filepath > $folder/$gsm_id.processed.txt
