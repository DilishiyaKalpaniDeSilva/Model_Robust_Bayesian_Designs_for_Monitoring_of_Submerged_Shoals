#!/bin/bash
path="codes/utility/main_utility_random_designs_GAMM.R"
log_path="sub_files/utility/log_files/exp_utility/"

B=125
B1=250

for indY in 1
do
  for indDC in {1..4}
  do
    for indD in {1..25}
    do
      one=$(qsub -N util_Y${indY}_DC${indDC}_D${indD} -v indexD=$indD,indexDC=$indDC,indexY=$indY,B=$B,B1=$B1,path=$path,log_path=$log_path sub_util_random_GAMM.sub)
      echo $one    
    done
  done
done






