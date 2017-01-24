#!/bin/bash

L="200"
q="1000"
temp="1.0"
lam="1.0"
b="1.0"
D="0.1"
t="300000"
teq="1000"
phi="1.0"
fracOccup="1.0"

a=$1
p=$2
run=$3
output_dir=$4
machine=$5

name="ratio_${L}_${L}_${q}_a_${a}_lam_${lam}_P_${p}_D_${D}_t_${t}_run_${run}"

log_file="${name}_${machine}.log"

nohup java cpm_boundary.CellPottsModel $L $L $q $temp $lam $a $b $p $phi $fracOccup $D $t $teq $run $output_dir &> $log_file &
 
