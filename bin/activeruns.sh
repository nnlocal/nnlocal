#!/bin/bash

OS=`uname -s`
if   [[ $OS == *"Linux"* ]]; then
    top -b n1 | grep nnlocal | wc -l
elif [[ $OS == *"Darwin"* ]]; then
    top -l 1 | grep nnlocal | wc -l
fi

exit

