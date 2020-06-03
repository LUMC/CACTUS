#!/bin/bash

declare -a arr=("S12118" "S144")

for i in "${arr[@]}"
do
   echo $i
   Rscript cactus_run.R  $1  $2  $3  $4  $i $5 &
done
wait
echo "All done"
