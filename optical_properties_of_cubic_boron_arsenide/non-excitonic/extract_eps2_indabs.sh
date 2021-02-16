#!/bin/bash

#1st column is energy 
#8th column corresponds to broadening of 0.1 eV, matching that used for direct
awk '{print $1,$8*4}' epsilon2_indabs.dat > indirect_eps2_col8.dat
#multiply by 4 corrects the bug in epw: Indabs.f90 factor of 0.5 in lines 282-288 is incorrect – holdover from unit conversions
