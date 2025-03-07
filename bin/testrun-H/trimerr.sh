#!/bin/bash
# script to reduce miss-binning using mergedata
# After ordering the entries in each bin by their error
# the first argument is the number of the entries to skip from above
# while the second count the entries to skip from below
(echo 1 $1 $2 ; ls -c1 ????-????.gnu ; echo "") | ../mergedata
rm fort.12
mv fort.13 nnlocal-$1-$2.top
exit

