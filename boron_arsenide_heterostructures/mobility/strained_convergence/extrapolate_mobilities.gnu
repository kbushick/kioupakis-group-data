file_name = 'convergence.dat'
xmin = 0.0
xmax_all = 2.5e-5
xmax_48 = 1.5e-5
xmax_72 = 3.5e-6

f(x) = a*x+b 
g(x) = c*x+d 
w(x) = m*x+n
y(x) = o*x+p

a = -4e7 
b = 2500
c = -4e7
d = 500
m = -4e7
n = 3600
o = -4e7 
p = 3100

fit [xmin:xmax_72] f(x) file_name using (1/$2):(($5+$6)/2) via a,b
fit [xmin:xmax_72] g(x) file_name using (1/$2):4 via c,d
fit [xmin:xmax_72] w(x) file_name using (1/$2):(($9+$10)/2) via m,n
fit [xmin:xmax_72] y(x) file_name using (1/$2):8 via o,p

print 'e mob low: ', d
print 'e mob high: ', b
print 'h mob low: ', p
print 'h mob high: ', n
plot file_name u (1/$2):4 w lp, file_name u (1/$2):(($5+$6)/2) w lp, f(x), g(x), file_name u (1/$2):8 w lp, file_name u (1/$2):(($9+$10)/2) w lp, w(x), y(x)
