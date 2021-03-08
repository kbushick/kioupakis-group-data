set style data dots
set nokey
set xrange [0: 4.41581]
set yrange [ -6.86312 : 43.30955]
set arrow from  1.16813,  -6.86312 to  1.16813,  43.30955 nohead
set arrow from  1.75219,  -6.86312 to  1.75219,  43.30955 nohead
set arrow from  2.16519,  -6.86312 to  2.16519,  43.30955 nohead
set arrow from  3.40418,  -6.86312 to  3.40418,  43.30955 nohead
set xtics ("G"  0.00000,"X"  1.16813,"W"  1.75219,"K"  2.16519,"G"  3.40418,"L"  4.41581)
 plot "Si_band.dat"
