#!/bin/bash

rm $3-plots.pdf
../genplots.sh $1 $2 $3
gnuplotsplit.sh genplots.gp
for file in *.eps
do
    epspdf $file
done
rm *.eps
pdftk *pdf cat output plots.pdf
rm $3-*.pdf
mv plots.pdf $3-plots.pdf
exit


