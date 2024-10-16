#!/bin/bash

# Usage: It expects a single file name or a list with usual * wildcard
# Run:
# > ../gengrids.sh *grid

for file in $@
do
    if [[ "$file" != "cfmc"*"xg"* ]];then
	../gridplot $file
	gnuplot $file.gnu
    fi
done

exit

