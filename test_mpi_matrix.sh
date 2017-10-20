#!/bin/bash
for ((n=6; n<=6;n++)) ; do for ((m=1;m<=$n;m+=1)) ; do for ((k=1;k<=$n;k++)) ; do echo "n=$n m=$m k=$k" ; mpirun -np $k ./a.out $n $m a.txt ; done ; done ; done
