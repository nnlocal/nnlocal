#!/bin/bash

> Timings.txt

seed=100
ncores=8
nprocessesgrid=8
nprocessesaccu=16
alpha=0.065
mtrim=$(echo "scale=4; $alpha * $nprocessesaccu" | bc)
mtrim=$(echo "scale=0; $mtrim / 1" | bc)

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
k=0
# two stages of importance sampling grid calculation
for igrid in `seq $maxgrid`
do

(echo -n st1 xg$igrid ' ' ; date ) >> Timings.txt
for j in `seq $nprocessesgrid`
do
    ((k++))
    ../nnlocal input.DAT $(( $seed + $k )) $j xg$igrid $nprocessesgrid $igrid > run-st1-xg$igrid-$j.log 2>&1 &
    limit_procs
done

wait

./reslist.sh st1-xg$igrid xg$igrid &

done

wait


(echo -n st2 '     ' ; date ) >> Timings.txt
for j in `seq $nprocessesaccu`
do
    ((k++))
    ../nnlocal input.DAT $(( $seed + $k )) $j st2 $nprocessesgrid $maxgrid > run-st2-$j.log 2>&1 &
    limit_procs
done
wait

./reslist.sh st2 st2 &


(echo 1 0 0 ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
mv fort.12 nnlocal-1.top

wait

(echo 2 0 0 ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
mv fort.12 nnlocal-2.top

wait

if [ $mtrim -eq 0 ]; then
    cp nnlocal-1.top nnlocal-0-0.top
else
    (echo 1 $mtrim $mtrim ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
    rm fort.12
    mv fort.13 nnlocal-$mtrim-$mtrim.top
fi

(echo -n end '     ' ; date ) >> Timings.txt
