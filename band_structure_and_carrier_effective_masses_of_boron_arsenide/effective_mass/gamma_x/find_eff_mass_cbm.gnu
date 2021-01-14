file_name = 'slice9_long.csv' 

lattice_param = 9.028354883
e_naught1 = 9.39222117602717 
xmin = -0.05
xmax = 0.05

x_shift = 1.02580480000000

inv_ang_to_inv_bohr = 0.529177249
ev_to_hartree = 1/27.21139614
qe_x_to_fit = 2*pi/lattice_param

f(x) = (-1+sqrt(1+4*a1*x**2/(2*m1)))/(2*a1)+e_naught1*ev_to_hartree

fit [xmin:xmax] f(x) file_name using ((($1-x_shift)*inv_ang_to_inv_bohr)):($2*ev_to_hartree) via a1,m1

plot [xmin:xmax] f(x), file_name using ((($1-x_shift)*inv_ang_to_inv_bohr)):($2*ev_to_hartree)
