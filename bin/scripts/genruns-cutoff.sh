#!/bin/bash

vals=( "1d-5" "2d-5" "5d-5" "1d-4" "2d-4" "5d-4" "1d-3" "2d-3" "5d-3" "1d-2" )

for i in `seq 0 9`
do
    echo $i
    for j in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
    do
	echo $j
	dir=testrun-cut-$i-seed-$j
	echo $dir
      	rm -fr $dir
	cp -r testrun-H-dummy $dir
	cd $dir
	ss=$((2089+10*$i+$j))
	echo $ss
	sed -i "s/seed/$ss/ ; s/valtiny/${vals[$i]}/" input.DAT
	../nnlocal > output.log &
	cd -
    done
done


exit

