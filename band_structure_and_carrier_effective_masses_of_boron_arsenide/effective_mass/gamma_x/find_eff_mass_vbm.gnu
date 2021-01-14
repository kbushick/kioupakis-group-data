file_name = 'slice_4.csv'
file_name2 = 'slice_6.csv'
file_name3 = 'slice_8.csv'

lattice_param = 9.028354883
e_naught4 = 7.13770754030170 
e_naught6 = 7.34345938062881
e_naught8 = 7.34345955276091
xmin = 0.0
xmax = 0.15

inv_ang_to_inv_bohr = 0.529177249
ev_to_hartree = 1/27.21139614
qe_x_to_fit = 2*pi/lattice_param

f(x) = (1-sqrt(1+4*a1*x**2/(2*m1)))/(2*a1)+e_naught4*ev_to_hartree
g(x) = (1-sqrt(1+4*a2*x**2/(2*m2)))/(2*a2)+e_naught6*ev_to_hartree
h(x) = (1-sqrt(1+4*a3*x**2/(2*m3)))/(2*a3)+e_naught8*ev_to_hartree


fit [xmin:xmax] f(x) file_name using ($1*inv_ang_to_inv_bohr):($2*ev_to_hartree) via a1,m1
fit [xmin:xmax] g(x) file_name2 using ($1*inv_ang_to_inv_bohr):($2*ev_to_hartree) via a2,m2
fit [xmin:xmax] h(x) file_name3 using ($1*inv_ang_to_inv_bohr):($2*ev_to_hartree) via a3,m3

plot [xmin:xmax] f(x), file_name using ($1*inv_ang_to_inv_bohr):($2*ev_to_hartree), g(x), file_name2 using ($1*inv_ang_to_inv_bohr):($2*ev_to_hartree), h(x), file_name3 using ($1*inv_ang_to_inv_bohr):($2*ev_to_hartree)

