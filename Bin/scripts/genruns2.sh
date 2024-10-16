#!/bin/bash

for j in `seq 100`
do
    dir=testrun-$j
    echo $dir
    rm -fr $dir
    cp -r testrun-H-dummy $dir
    cd $dir
    ss=$((1089+$j))
    sed -i "s/seed/$ss/" input.DAT
    ../nnlocal2 > output.log &
    cd -
done

exit

