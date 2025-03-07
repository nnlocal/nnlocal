#!/bin/bash

> Timings.txt

ncores=8
nprocessesgrid=8
nprocessesaccu=128

# First compile the nnlocal executable in the ../../ directory
#

# the following function limits the number of subprocesses
# to be not larger than the number of cores specified by the
# user

function limit_procs {
    while [ `jobs -p | wc -w` -ge $ncores ]
    do
	sleep 1
    done
}

maxgrid=4
# two stages of importance sampling grid calculation
for igrid in `seq $maxgrid`
do

(echo -n st1 xg$igrid ' ' ; date ) >> Timings.txt

for i in `seq $nprocessesgrid`
do
    ../nnlocal input.DAT $i xg$igrid- $nprocessesgrid $igrid > run-st1-xg$igrid-$i.log 2>&1 &
    limit_procs
done

wait

./reslist.sh st1-xg$igrid xg$igrid &

done

wait


(echo -n st2 '     ' ; date ) >> Timings.txt
for i in `seq $nprocessesaccu`
do
    ../nnlocal input.DAT $i st2- $nprocessesgrid $maxgrid > run-st2-$i.log 2>&1 &
    limit_procs
done
wait

./reslist.sh st2 st2 &


(echo 1 0 0 ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
mv fort.12 nnlocal-1.top

wait

(echo 2 0 0 ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
mv fort.12 nnlocal-2.top

(echo -n end '     ' ; date ) >> Timings.txt

