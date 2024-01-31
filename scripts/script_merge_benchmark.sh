#!/bin/bash/env sh
remove=$(realpath $1)
all=$(realpath $2)
output=$(realpath $3)

#remove="/net/beegfs/cfg/tgac/ferdinand/WES-snake/benchmarks/104-026-9_HGTGWDSXY_S33/"
#names=$(ls /net/beegfs/cfg/tgac/ferdinand/WES-snake/benchmarks/104-026-9_HGTGWDSXY_S33/*)

names=$(ls $all )
one_file=$(echo $names | awk '{print $1}')
header=$(cat ${one_file} |  egrep  'h:m:s')

printf "RULE\t%s\n" "$header" > ${output} #docu.tsv
for file in ${names[@]}; 
do 
   name=${file/#$remove"/"}   # ook normaal
   name2=$(echo $name | sed 's/.bwa.benchmark.txt//g')   #ook normaal
   data=$(cat $file | egrep -v 'h:m:s')   #dit is raarder
   printf "%s\t%s\n" "$name2" "$data" >> ${output} #docu.tsv 
done


