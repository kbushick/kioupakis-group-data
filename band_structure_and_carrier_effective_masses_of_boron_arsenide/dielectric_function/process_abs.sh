#!/bin/bash
num=18
tail -2601 absorption_noeh_"$num".dat | awk '{print $1, $2}' > temp
tail -2601 absorption_eh_"$num".dat | awk '{print $2}' > temp1
paste temp temp1 | awk '{print $0}' > abs_"$num".dat
rm temp temp1
