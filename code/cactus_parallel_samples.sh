#!/bin/bash

declare -a arr=("S12118" "S144")

for i in "${arr[@]}"
do
   echo $i
   Rscript cactus_run.R  "../input/WES/wes.rds"  "../input/scRNA/ac.rds"  "../results/"   "../tree/"  $i "../input/scBCR_GEX.rds"  "0.5" "9.5" &
done
wait
echo "All done"
