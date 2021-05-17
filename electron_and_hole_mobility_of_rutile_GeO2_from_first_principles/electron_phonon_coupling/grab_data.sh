#!/bin/bash

#grep "iq = " epw2.out -A 21 > tmp
grep "iq = " epw2.out -A 2181 > tmp # modes*bnds*bnds+3 (18*(12-4+1)*(12-4+1))+3
echo "--" >> tmp

