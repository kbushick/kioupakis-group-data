#!/bin/bash

for ctype in lda
do

#valence bands
head -643  pz_$ctype.dat > lda/band1_$ctype.csv
head -643  pz_so_$ctype.dat > lda_so/band1_so_$ctype.csv
head -1286 pz_$ctype.dat | tail -643 > lda/band2_$ctype.csv
head -1286 pz_so_$ctype.dat | tail -643 > lda_so/band2_so_$ctype.csv
head -1929 pz_$ctype.dat | tail -643 > lda/band3_$ctype.csv
head -1929 pz_so_$ctype.dat | tail -643 > lda_so/band3_so_$ctype.csv
head -2572 pz_$ctype.dat | tail -643 > lda/band4_$ctype.csv
head -2572 pz_so_$ctype.dat | tail -643 > lda_so/band4_so_$ctype.csv
head -3215 pz_so_$ctype.dat | tail -643 > lda_so/band5_so_$ctype.csv
head -3858 pz_so_$ctype.dat | tail -643 > lda_so/band6_so_$ctype.csv
head -4501 pz_so_$ctype.dat | tail -643 > lda_so/band7_so_$ctype.csv
head -5144 pz_so_$ctype.dat | tail -643 > lda_so/band8_so_$ctype.csv

#conduction bands
head -3215 pz_$ctype.dat | tail -643 > lda/band5_$ctype.csv
head -5787 pz_so_$ctype.dat | tail -643 > lda_so/band9_so_$ctype.csv
head -3858 pz_$ctype.dat | tail -643 > lda/band6_$ctype.csv
head -6430 pz_so_$ctype.dat | tail -643 > lda_so/band10_so_$ctype.csv
head -4501 pz_$ctype.dat | tail -643 > lda/band7_$ctype.csv
head -7073 pz_so_$ctype.dat | tail -643 > lda_so/band11_so_$ctype.csv
head -5144 pz_$ctype.dat | tail -643 > lda/band8_$ctype.csv
head -7716 pz_so_$ctype.dat | tail -643 > lda_so/band12_so_$ctype.csv
head -8359 pz_so_$ctype.dat | tail -643 > lda_so/band13_so_$ctype.csv
head -9002 pz_so_$ctype.dat | tail -643 > lda_so/band14_so_$ctype.csv
head -9645 pz_so_$ctype.dat | tail -643 > lda_so/band15_so_$ctype.csv
head -10288 pz_so_$ctype.dat | tail -643 > lda_so/band16_so_$ctype.csv

done

for ctype in gw
do

#valence bands
head -643  pz_$ctype.dat > gw/band1_$ctype.csv
head -643  pz_so_$ctype.dat > gw_so/band1_so_$ctype.csv
head -1286 pz_so_$ctype.dat | tail -643 > gw_so/band2_so_$ctype.csv
head -1929 pz_$ctype.dat | tail -643 > gw/band2_$ctype.csv
head -1929 pz_so_$ctype.dat | tail -643 > gw_so/band3_so_$ctype.csv
head -2572 pz_so_$ctype.dat | tail -643 > gw_so/band4_so_$ctype.csv
head -3215 pz_$ctype.dat | tail -643 > gw/band3_$ctype.csv
head -3215 pz_so_$ctype.dat | tail -643 > gw_so/band5_so_$ctype.csv
head -3858 pz_so_$ctype.dat | tail -643 > gw_so/band6_so_$ctype.csv
head -4501 pz_$ctype.dat | tail -643 > gw/band4_$ctype.csv
head -4501 pz_so_$ctype.dat | tail -643 > gw_so/band7_so_$ctype.csv
head -5144 pz_so_$ctype.dat | tail -643 > gw_so/band8_so_$ctype.csv

#conduction bands
head -5787 pz_$ctype.dat | tail -643 > gw/band5_$ctype.csv
head -5787 pz_so_$ctype.dat | tail -643 > gw_so/band9_so_$ctype.csv
head -6430 pz_so_$ctype.dat | tail -643 > gw_so/band10_so_$ctype.csv
head -7073 pz_$ctype.dat | tail -643 > gw/band6_$ctype.csv
head -7073 pz_so_$ctype.dat | tail -643 > gw_so/band11_so_$ctype.csv
head -7716 pz_so_$ctype.dat | tail -643 > gw_so/band12_so_$ctype.csv
head -8359 pz_$ctype.dat | tail -643 > gw/band7_$ctype.csv
head -8359 pz_so_$ctype.dat | tail -643 > gw_so/band13_so_$ctype.csv
head -9002 pz_so_$ctype.dat | tail -643 > gw_so/band14_so_$ctype.csv
head -9645 pz_$ctype.dat | tail -643 > gw/band8_$ctype.csv
head -9645 pz_so_$ctype.dat | tail -643 > gw_so/band15_so_$ctype.csv
head -10288 pz_so_$ctype.dat | tail -643 > gw_so/band16_so_$ctype.csv

done
