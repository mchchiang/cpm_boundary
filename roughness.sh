#!/bin/bash

temp="1.0"
lam="1.0"
a="2.0"
b="2.0"
p="10.0"
D="0.1"
t="1000500"
teq="500"
fracOccup="0.5"

L=$1
q=$2
phi=$3

run=$4
output_dir=$5
machine=$6

name="rough_${L}_${L}_${q}_a_${a}_b_${b}_lam_${lam}_P_${p}_D_${D}_fracMo_${phi}_fracOccup_${fracOccup}_t_${t}_run_${run}"

log_file="${name}_${machine}.log"

nohup java cpm_boundary.CellPottsModel $L $L $q $temp $lam $a $b $p $phi $fracOccup $D $t $teq $run $output_dir &> $log_file &
 
