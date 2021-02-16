awk '{print $5}' avg_LS.dat > temp_eps2_dir.dat
paste indirect_eps2_col8_interpolated.dat temp_eps2_dir.dat | awk '{print $0, $2+$3}' > all_eps2.dat
awk '{print $1,$4}' all_eps2.dat > total_eps2.dat  
./kk > total_eps1.dat
awk '{print $2}' total_eps2.dat > temp_eps2_tot_col
#put x (energy) eps1 eps2 in same file
paste total_eps1.dat temp_eps2_tot_col | column -t > temp_eps1_eps2
#calculate n and k using equations 1.22, 1.23 from "Optical Properties of Solids" 
awk '{print $0, 1 / sqrt(2) * sqrt($2 + sqrt($2 * $2 + $3 * $3)), 1 / sqrt(2) * sqrt(-1 * $2 + sqrt($2 * $2 + $3 * $3))}' temp_eps1_eps2 > temp_1_2_n_k 
#calcuate absorption using eq 1.16, in units of inv meters
awk '{print $0, 4*3.14159265 * $5 * $1 / 1.239841973862093e-06}' temp_1_2_n_k | column -t > optical_data_total.dat
#columns have energy, eps1, eps2, n, k, alpha in columns
rm temp*
