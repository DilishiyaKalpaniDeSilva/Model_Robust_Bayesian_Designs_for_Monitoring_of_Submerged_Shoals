#!/bin/bash
path1="codes/utility/main_relative_utility_pol_without_yre.R"
path2="codes/utility/main_relative_utility_pol_with_yre.R"
log_path="sub_files/utility/log_files/rel_utility/"

B=2
B1=50

for indY in {1..3}
do
  for indO in {1..3}
  do
    one=$(qsub -N util_Y${indY}_O${indO} -v indexO=$indO,indexY=$indY,B=$B,B1=$B1,path=$path1,log_path=$log_path sub_rel_utility.sub)
    echo $one    
  done
done

for indY in 4
do
  for indO in {1..3}
  do
    two=$(qsub -N util_Y${indY}_O${indO} -v indexO=$indO,indexY=$indY,B=$B,B1=$B1,path=$path2,log_path=$log_path sub_rel_utility.sub)
    echo $two    
  done
done




