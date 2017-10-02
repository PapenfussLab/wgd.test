#!/bin/bash

NS="0 1 2 3 4 5"
MS="-1 0 1 2 3"
BOOLS="True False"
B="True"

for N in $NS;
do
  for M in $MS;
  do
    for FA in $BOOLS; 
    do
      echo $N $M $B $FA
      /Users/lmcintosh/Documents/SageMath/sage GD_v1.0.sage $N $M $B $FA
    done;
  done;
done;

#for N in $NS;
#do
#  for M in $MS;
#  do
#    for B in $BOOLS;
#    do
#      for FA in $BOOLS; 
#      do
#        echo $N $M $B $FA
#        /Users/lmcintosh/Documents/SageMath/sage GD_v1.0.sage $N $M $B $FA
#      done;
#    done;
#  done;
#done;
