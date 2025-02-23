#!/bin/bash
path="codes/models/main_model_selection.R"
log_path="sub_files/models/log_files/model_selection/"

from=1
to=241	# no of models

for indY in {1..4} 	# index of the data set
do
  job_u=$(qsub -N Model_selection${indY} -v indexY=$indY,path=$path,log_path=$log_path -J $from-$to sub_model_select.pbs)
  echo $job_u
done
