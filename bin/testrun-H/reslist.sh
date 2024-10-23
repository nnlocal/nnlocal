#!/bin/bash

nfiles=`ls -c1 run-$1-*.log | wc -l`

outfile=results-$1.dat

echo "$2 = {" > $outfile  

for file in run-$1-*.log
do    
    read -ra ADDR <<< "$(grep Value $file)"
    if [ ! -z "$ADDR" ]
    then
	echo $file
        echo '{' ${ADDR[6]} ',' ${ADDR[8]} '},' >> $outfile
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

exit

