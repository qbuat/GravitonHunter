#!/bin/bash

#for file in `ls ~/mysps/DATASETS/data12_8TeV.00201190*debug*/*`
for file in `ls ~/mysps/DATASETS/data12_8TeV.*debug*/*`
do
  name=`echo $file`
  run=`echo $name | awk -F '.' '{print $3}' | sed -r 's/^.{2}//'`
  echo $file >> input_debug_run${run}.txt
done