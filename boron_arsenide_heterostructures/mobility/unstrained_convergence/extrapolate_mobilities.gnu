file_name = 'convergence.dat'

xmin = 0.0
xmax_all = 2.5e-5
xmax_48 = 1.5e-5
xmax_72 = 3.5e-6

f(x) = a*x+b 
g(x) = c*x+d 

a = 6e7
b = 1250
c = -4e7
d = 1390

fit [xmin:xmax_72] f(x) file_name using (1/$2):7 via a,b
fit [xmin:xmax_72] g(x) file_name using (1/$2):11 via c,d

print 'e mob avg: ', b
print 'h mob avg: ', d
plot file_name u (1/$2):7 w lp, f(x), file_name u (1/$2):11 w lp, g(x)
