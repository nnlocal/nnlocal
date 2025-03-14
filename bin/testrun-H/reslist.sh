#!/bin/bash

nfiles=`ls -c1 run-$1-*.log | wc -l`

outfile=results-$1.dat

val=0
err=0
scale=6

echo "$2 = {" > $outfile  

for file in run-$1-*.log
do    
    read -ra ADDR <<< "$(grep Value $file)"
    if [ ! -z "$ADDR" ]
    then
	echo $file
        echo '{' ${ADDR[6]} ',' ${ADDR[8]} '},' >> $outfile
	val=$(echo "scale=$scale; $val + ${ADDR[6]}" | bc)
	err=$(echo "scale=$scale; $err + ${ADDR[8]}^2" | bc)
#	echo $val $err
    fi
done

OS=`uname -s`
if   [[ $OS == *"Linux"* ]]; then
    sed -i '$s/}\,/}};/g' $outfile
elif [[ $OS == *"Darwin"* ]]; then
    sed -i '' '$s/}\,/}};/g' $outfile
fi

data=`cat $outfile | wc -l`

echo 'Collected' $(( data - 1 )) '/' $nfiles 'results'
val=$(echo "scale=$scale; $val/($data-1)" | bc)
err=$(echo "scale=$scale; sqrt($err)/($data-1)" | bc)
echo 'Value of cross section is: ' $val ' +/- ' $err 'fb'
exit

