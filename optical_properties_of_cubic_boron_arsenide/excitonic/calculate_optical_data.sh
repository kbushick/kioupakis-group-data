rm eps2.dat
awk '{print $2}' absorption_eh.dat > eps2col.tmp
awk '{print $1,$2}' absorption_eh.dat > eps2.dat 
./kk > "eh"_eps1.dat
#put x (energy) eps1 eps2 in same file
paste "eh"_eps1.dat eps2col.tmp | column -t > eps1_eps2.tmp
#calculate n and k using equations 1.22, 1.23 from "Optical Properties of Solids" 
awk '{print $0, 1 / sqrt(2) * sqrt($2 + sqrt($2 * $2 + $3 * $3)), 1 / sqrt(2) * sqrt(-1 * $2 + sqrt($2 * $2 + $3 * $3))}' eps1_eps2.tmp > e_1_2_n_k.tmp 
#calcuate absorption using eq 1.16, in units of inv meters
awk '{print $0, 4*3.14159265 * $5 * $1 / 1.239841973862093e-06}' e_1_2_n_k.tmp | column -t > optical_data_"eh".dat
#columns have energy, eps1, eps2, n, k, alpha in columns
rm *tmp

