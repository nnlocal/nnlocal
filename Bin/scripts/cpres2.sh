#!/bin/bash


mkdir newresults

for i in `seq 100`
do
    echo $i
    dir=testrun-$i
    mkdir newresults/$dir
    cp $dir/output.log newresults/$dir/nohup.out
    done
done

cp genruns2.sh cpres2.sh newresults/

exit

