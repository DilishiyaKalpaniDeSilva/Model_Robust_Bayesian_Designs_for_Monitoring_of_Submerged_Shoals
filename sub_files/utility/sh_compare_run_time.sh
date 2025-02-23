#!/bin/bash
path="codes/utility/main_utility_run_time_GAMM.R"
log_path="sub_files/utility/log_files/run_time/"

from=1
to=9

for indY in 1
do 
  for indD in {1..50}
  do
    for indB in {1..2}
    do
      for indB1 in {1..2}
      do
		    one=$(qsub -N run_time_Y${indY}_D${indD}_B${indB}_B1${indB1} -v indexY=$indY,indexD=$indD,indexB=$indB,indexB1=$indB1,path=$path,log_path=$log_path -J $from-$to sub_compare_run_time.sub)
		    echo $one 
      done
    done
  done
done


