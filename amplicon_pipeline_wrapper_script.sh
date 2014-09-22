#!/bin/bash

fd=$1

read -ra r1 <<< `find $fd | grep "R1"`
read -ra r2 <<< `find $fd | grep "R2"`

for ((i=0;i<${#r1[@]};++i)); do
    #hacky way of getting kk # name
    name=`echo ${r1[i]} | grep -o '[^/]*$' | cut -d "_" -f 1 | cut -d "-" -f 2`
    python amplicon_pipeline.py --name $name --fwd ${r1[i]} --rev ${r2[i]}
done
