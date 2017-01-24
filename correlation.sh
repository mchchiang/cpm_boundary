#!/bin/bash

temp="1.0"
lam="1.0"
a="2.0"
b="2.0"
p="10.0"
D="0.1"
t="1000500"
teq="500"
fracOccup="0.50"

L=$1
q=$2
phi=$3

run=$4
output_dir=$5
machine=$6

name="rough_${L}_${L}_${q}_a_${a}_b_${b}_lam_${lam}_P_${p}_D_${D}_fracMo_${phi}_fracOccup_${fracOccup}_t_${t}_run_${run}"

log_file="${name}_correlation_${machine}.log"

nohup java cpm_boundary.Correlation 3 $t 0 100000 "${output_dir}/${name}.dat" "${output_dir}/${name}_correlation.dat" &> $log_file &
