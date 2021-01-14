for band in 4 6 8
do
head -20 band"$band"_so_gw.csv > gamma_x/slice_$band.csv
head -311 band"$band"_so_gw.csv | tail -20 > gamma_l/slice_$band.csv
head -292 band"$band"_so_gw.csv | tail -20 > gamma_k/slice_$band.csv
done

for band in 3 5 7
do
head -292 band"$band"_so_gw.csv | tail -20 > gamma_k/slice_$band.csv
done
