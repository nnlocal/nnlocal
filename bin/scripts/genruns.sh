#!/bin/bash

vals1=( "125d0" "400d0" "1000d0")

vals2=( "10d0" "20d0" )

for i in `seq 0 2`
do
    for k in `seq 0 1`
    do
	for j in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	do
	    echo $j
	    dir=testrun-mass-$i-cut-$k-seed-$j
	    echo $dir
      	    rm -fr $dir
	    cp -r testrun-H-dummy $dir
	    cd $dir
	    ss=$((1089+10*$i+$j))
	    sed -i "s/seed/$ss/ ; s/mmm/${vals1[$i]}/ ; s/ppp/${vals2[$k]}/" input.DAT
	    ../nnlocal > output.log &
	    cd -
	done
    done
done

exit

