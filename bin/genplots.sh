#!/bin/bash
if ! [ -e "$1" ]
then
echo "usage: genplots.sh file1 file2 nameoutput"
exit -1
fi
if ! [ -e "$2" ]
then
echo "usage: genplots.sh file1 file2 nameoutput"
exit -1
fi
if [ z$3 = z ]
then
echo "usage: genplots.sh file1 file2 nameoutput"
exit -1
fi
>genplots.gp
i=0
fil1=$1
fil2=$2
cat $fil2 | grep '#' | while read line
do
title=`echo "'""$line""'"|sed 's/#//'`
echo "$title"

ttt=`echo "$line"|sed 's/#//g ; s/ /_/g ; s/\//o/g ' `

file=$3-$ttt

cat <<EOF >> genplots.gp
reset
set term post eps enhanced color dashed "Helvetica" 24
set output '$file.eps'
iindex = $i
mytitle = $title
infil1 = "$fil1"
infil2 = "$fil2"
both = '< paste $fil1 $fil2'
set multiplot

set origin 0.0,0.3
set size 1.0,0.7
set lmargin at screen 0.2
set rmargin at screen 0.95
set tmargin at screen 0.90
set bmargin at screen 0.4

#set log x
set log y

set datafile fortran
set style data xyerrorlines

set xtics format ""


EOF
cat <<'EOF' >> genplots.gp
set title mytitle
set xrange [] writeback
plot infil1 index  iindex using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title infil1,\
    infil2 index  iindex  using ($1+$2)/2:3:1:2:($3-$4):($3+$4) title infil2

set tmargin at screen 0.4
set bmargin at screen 0.1

set nolog y
set yrange [0:2]
set xrange restore
set format x
unset title
plot both index iindex \
using ($1+$2)/2:($3/$7):1:2:($3/$7)*(1-(($4/$3)**2+($8/$7)**2)**0.5):($3/$7)*(1+(($4/$3)**2+($8/$7)**2)**0.5) title "Ratio", 1 with lines linestyle -1

unset multiplot

EOF

i=$[$i+1]
done
