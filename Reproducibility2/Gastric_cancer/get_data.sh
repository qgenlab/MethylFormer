


path_to_fasterq_dump=$1

while read srr; do [ -z "$srr" ] && continue; echo "Downloading $srr"; $path_to_fasterq_dump/fasterq-dump --split-files $srr; done < SRR_Acc_List.txt


