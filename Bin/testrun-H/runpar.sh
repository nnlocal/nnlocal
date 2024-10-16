#!/bin/bash

> Timings.txt

ncores=100
nprocessesgrid=100
nprocessesaccu=400
# First compile the pwhg_main executable in the ../ directory
#

# the following function limits the number of subprocesses
# to be not larger than the number of cores specified by the
# user

function limit_procs {
    while [ `jobs -p | wc -w` -gt $ncores ]
    do
	sleep 1
    done
}

maxgrid=5
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

done

wait


(echo -n st2 '     ' ; date ) >> Timings.txt
for i in `seq $nprocessesaccu`
do
    ../nnlocal input.DAT $i st2- $nprocessesgrid $maxgrid > run-st2-$i.log 2>&1 &
    limit_procs
done
wait

(echo 1 ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
mv fort.12 nnlocal-1.top

wait

(echo 2 ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
mv fort.12 nnlocal-2.top

(echo -n end '     ' ; date ) >> Timings.txt

