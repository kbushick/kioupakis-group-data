set style data dots
set nokey
set xrange [0: 4.41581]
set yrange [ -7.62313 : 41.88412]
set arrow from  1.16813,  -7.62313 to  1.16813,  41.88412 nohead
set arrow from  1.75219,  -7.62313 to  1.75219,  41.88412 nohead
set arrow from  2.16519,  -7.62313 to  2.16519,  41.88412 nohead
set arrow from  3.40418,  -7.62313 to  3.40418,  41.88412 nohead
set xtics ("G"  0.00000,"X"  1.16813,"W"  1.75219,"K"  2.16519,"G"  3.40418,"L"  4.41581)
 plot "Si_band.dat"
