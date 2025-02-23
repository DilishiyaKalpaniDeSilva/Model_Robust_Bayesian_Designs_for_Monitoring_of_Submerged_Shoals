#!/bin/bash
path1="codes/models/main_priors_GAMM.R"
path2="codes/models/main_priors_polynomial.R"

log_path="sub_files/models/log_files/models/priors/"

indF=1	# index of the fishnet size
from=1
to=4

## run GAMM models
job_u=$(qsub -N priors_GAMM_F${indF} -v indexF=$indF,path=$path1,log_path=$log_path -J $from-$to sub_priors_GAMM.sub)
echo $job_u

## run polynomial models
for indO in {1..3}  # index for the degree of polynomial
do
  job_u=$(qsub -N priors_pol_F${indF}_O${indO} -v indexF=$indF,indexO=$indO,path=$path2,log_path=$log_path -J $from-$to sub_priors_pol.sub)
  echo $job_u
done

