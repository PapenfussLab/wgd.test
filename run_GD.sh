#!/bin/sh

minimum_n=0
maximum_n=3
#maximum_n=5
minimum_m=-1
maximum_m=3
#maximum_m=4

for n in $(seq "$minimum_n" "$maximum_n"); do
  for m in $(seq "$minimum_m" "$maximum_m"); do
      printf 'N=%d; M=%d\n' "$n" "$m"
      sage GD_v1.0.sage --balanced $n $m
      sage GD_v1.0.sage --balanced --forced_alive $n $m
  done;
done;
