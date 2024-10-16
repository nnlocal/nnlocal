#!/bin/bash

vals=( "1d-5" "2d-5" "5d-5" "1d-4" "2d-4" "5d-4" "1d-3" "2d-3" "5d-3" "1d-2" )

mkdir results

for i in `seq 0 9`
do
    echo $i
    for j in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
    do
	echo $j
	dir=testrun-cut-$i-seed-$j
	echo $dir
      	mkdir results/$dir
	cp $dir/nohup.out results/$dir/nohup.out
    done
done

cp genruns.sh cpres.sh results/

exit

