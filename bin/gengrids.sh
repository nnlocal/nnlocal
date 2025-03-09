#!/bin/bash

# Usage: It expects a single file name or a list with usual * wildcard
# Run:
# > ../gengrids.sh *grid

for file in $@
do
    ../gridplot $file
    gnuplot $file.gnu
done

exit

