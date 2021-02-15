#!/bin/bash

rm convergence.dat

for f in 36 48 60 72 96 120 
do

f3=$(echo $f*$f*$f | bc)
irq=$(grep only $f/epw2.out | awk '{printf "%.4f", $6}')

ex=$(tail -1 $f/BAs_elcond_e | awk '{printf "%.4f", $3}')
ey=$(tail -1 $f/BAs_elcond_e | awk '{printf "%.4f", $7}')
ez=$(tail -1 $f/BAs_elcond_e | awk '{printf "%.4f", $11}')

hx=$(tail -1 $f/BAs_elcond_h | awk '{printf "%.4f", $3}')
hy=$(tail -1 $f/BAs_elcond_h | awk '{printf "%.4f", $7}')
hz=$(tail -1 $f/BAs_elcond_h | awk '{printf "%.4f", $11}')

emx=$(grep "Elec mob" $f/epw2.out -A 6 | tail -4 | grep x-axis | awk '{printf "%.4f", $4}')
emy=$(grep "Elec mob" $f/epw2.out -A 6 | tail -4 | grep y-axis | awk '{printf "%.4f", $1}')
emz=$(grep "Elec mob" $f/epw2.out -A 6 | tail -4 | grep z-axis | awk '{printf "%.4f", $1}')
ema=$(grep "Elec mob" $f/epw2.out -A 6 | tail -4 | grep avg | awk '{printf "%.4f", $1}')

hmx=$(grep "Hole mob" $f/epw2.out -A 6 | tail -4 | grep x-axis | awk '{printf "%.4f", $4}')
hmy=$(grep "Hole mob" $f/epw2.out -A 6 | tail -4 | grep y-axis | awk '{printf "%.4f", $1}')
hmz=$(grep "Hole mob" $f/epw2.out -A 6 | tail -4 | grep z-axis | awk '{printf "%.4f", $1}')
hma=$(grep "Hole mob" $f/epw2.out -A 6 | tail -4 | grep avg | awk '{printf "%.4f", $1}')

esprd=$(echo $emx - $emz | bc)
hsprd=$(echo $hmx - $hmz | bc)

echo $f $f3 $irq $emx $emy $emz $ema $hmx $hmy $hmz $hma $ex $ey $ez $hx $hy $hz $esprd $hsprd >> convergence.dat

done
