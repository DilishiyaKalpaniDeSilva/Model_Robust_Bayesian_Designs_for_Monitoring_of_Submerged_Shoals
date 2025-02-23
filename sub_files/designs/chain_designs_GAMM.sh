#!/bin/bash

# Code paths and log directory
#dpath="GAMM/without_yre"    #uncomment this to run GAMM designs without year random effect
dpath="GAMM/with_yre"     #uncomment this to run GAMM designs with year random effect
path0="codes/designs/${dpath}/0_initial_utility.R"
path1="codes/designs/${dpath}/1_approximate_utility_angle.R"
path2="codes/designs/${dpath}/2_remaining_angles.R"
path3="codes/designs/${dpath}/3_approximate_utility_angle_rest.R"
path4="codes/designs/${dpath}/4_optimize_angle_ACE.R"
path5="codes/designs/${dpath}/5_summarize_transect.R"
path7="codes/designs/${dpath}/7_compare_iteration.R"
log_path1="sub_files/designs/log_files/${dpath}/"

# Input parameters
B=250                 # no of Monte Carlo Draws for expected utility 
B1=500               # no of Monte Carlo Draws for Laplace approximation
from1=1; to1=20; from2=1; to2=10   # job array index for 20 angles
ncl=12               # no of clusters to use in foreach
indI=1              # initial design index

# Script names
sub_init="sub_init_util.pbs"
sub_angle="sub_util_angle.pbs"
sub_rest="sub_angle_rest.pbs"
sub_opt="sub_optimize_angle.pbs"
sub_sum="sub_sum_trans.pbs"
sub_com="sub_com_iter.pbs"


# Loop over indices and iterations
for indY in 1; do     # index for prior data set
  log_path="${log_path1}/${indY}/${indI}/"
  
  #Initial job submission for initialization
  #init=$(qsub -N init_util_Y${indY} -v log_path=$log_path,ncl=$ncl,indexY=$indY,indexI=$indI,B=$B,B1=$B1,path=$path0 $sub_init)
  #echo "Initial job submitted with ID: $init"

  for iter in 1; do     # loop over itertaions
    for indT in 18; do      # loop over transects
        
      jname_t="Y${indY}_I${indI}_T${indT}_Iter${iter}"
      common_params="log_path=$log_path,ncl=$ncl,indexI=$indI,indexY=$indY,iter=$iter,indexTr=$indT,B=$B,B1=$B1"
      
      if [ "$iter" -eq 1 ] && [ "$indT" -eq 1 ]; then
          dep_cond="-W depend=afterok:$init"  # First iteration of the first transect depends on the initial job
      elif [ "$iter" -eq 1 ] && [ "$indT" -eq 18 ]; then
          dep_cond=""  # No initial dependency for subsequent transects within the first iteration
      elif [ "$iter" -gt 1 ] && [ "$indT" -eq 1 ]; then
          dep_cond="-W depend=afterany:$one_c"  # No initial dependency for subsequent transects within the first iteration
      else
          dep_cond="-W depend=afterany:$one_t"  # Subsequent iterations depend on the completion of the composite job from the previous iteration
      fi

      dep_cond_t="depend=afterany"
      for indP in {1..22}; do     # loop over different angles
          jname_p="${jname_t}_P${indP}"
          params="${common_params},indexP=$indP"

          onea=$(qsub -N "util_angle_${jname_p}" -v $params,path=$path1 ${dep_cond} -J $from1-$to1 $sub_angle)
          echo "Submitted: $onea, Dependency: ${dep_cond}"  

          oner1=$(qsub -N "angle_rest_${jname_p}" -v $params,path=$path2 -W depend=afterany:$onea $sub_rest)
          echo "Submitted: $oner1, Dependency: $onea"  

          onea1=$(qsub -N "util_angle_rest_${jname_p}" -v $params,path=$path3 -W depend=afterany:$oner1 -J $from2-$to2 $sub_angle)
          echo "Submitted: $onea1, Dependency: $oner1"  

          onep=$(qsub -N "optimize_angle_${jname_p}" -v $params,path=$path4 -W depend=afterany:$onea1 $sub_opt)
          echo "Submitted: $onep, Dependency: $onea1"  
          
          dep_cond_t+=":$onep"
      done

      # Summarize at the end of each iteration for the transect
      one_t=$(qsub -N "sum_Trans_${jname_t}" -v $common_params,path=$path5 -W "${dep_cond_t}" $sub_sum)
      echo "Submitted: $one_t, Dependency: $dep_cond_t"  
      
    done
    
    one_c=$(qsub -N "com_iter_${jname_p}" -v $common_params,path=$path7 -W depend=afterany:$one_t $sub_com)
	echo $one_c
    
  done
done
