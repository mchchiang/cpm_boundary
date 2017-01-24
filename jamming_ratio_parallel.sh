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

start_a=$1
max_a=$2
inc_a=$3

start_p=$4
max_p=$5
inc_p=$6

trials=$7
cores=$8
output_dir=$9
machine=${10}

name="${L}_${L}_${q}_a_${start_a}-${max_a}-${inc_a}_lam_${lam}_P_${start_p}-${max_p}-${inc_p}_D_${D}_t_${t}_trials_${trials}"

log_file="ratio_${name}_${machine}.log"

nohup java cpm_boundary.CPMParallelRun $L $L $q $temp $lam $start_a $max_a $inc_a $b $start_p $max_p $inc_p $D $t $teq $trials $cores $output_dir &> $log_file &
 
