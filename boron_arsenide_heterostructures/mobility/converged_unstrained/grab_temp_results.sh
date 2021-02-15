#!/bin/bash

rm temp_72_data.*

for t in $(seq 100 50 500)
do

emf=$(grep final epw2.out -A 82 | tail -81 | grep Elec -A 40 | grep $t.000 -A 3 | grep x-axis | awk '{printf "%.4f", $2}')
emx=$(grep final epw2.out -A 82 | tail -81 | grep Elec -A 40 | grep $t.000 -A 3 | grep x-axis | awk '{printf "%.4f", $4}')
emy=$(grep final epw2.out -A 82 | tail -81 | grep Elec -A 40 | grep $t.000 -A 3 | grep y-axis | awk '{printf "%.4f", $1}')
emz=$(grep final epw2.out -A 82 | tail -81 | grep Elec -A 40 | grep $t.000 -A 3 | grep z-axis | awk '{printf "%.4f", $1}')
ema=$(grep final epw2.out -A 82 | tail -81 | grep Elec -A 40 | grep $t.000 -A 3 | grep avg | awk '{printf "%.4f", $1}')

hmf=$(grep final epw2.out -A 82 | tail -81 | grep Hole -A 40 | grep $t.000 -A 3 | grep x-axis | awk '{printf "%.4f", $2}')
hmx=$(grep final epw2.out -A 82 | tail -81 | grep Hole -A 40 | grep $t.000 -A 3 | grep x-axis | awk '{printf "%.4f", $4}')
hmy=$(grep final epw2.out -A 82 | tail -81 | grep Hole -A 40 | grep $t.000 -A 3 | grep y-axis | awk '{printf "%.4f", $1}')
hmz=$(grep final epw2.out -A 82 | tail -81 | grep Hole -A 40 | grep $t.000 -A 3 | grep z-axis | awk '{printf "%.4f", $1}')
hma=$(grep final epw2.out -A 82 | tail -81 | grep Hole -A 40 | grep $t.000 -A 3 | grep avg | awk '{printf "%.4f", $1}')


esprd=$(echo $emx - $emz | bc)
hsprd=$(echo $hmx - $hmz | bc)

echo $t $emx $emy $emz $ema $hmx $hmy $hmz $hma $emf $hmf $esprd $hsprd >> temp_72_data.dat
echo $t,$emx,$emy,$emz,$ema,$hmx,$hmy,$hmz,$hma,$emf,$hmf,$esprd,$hsprd >> temp_72_data.csv

done
