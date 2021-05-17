file_elec= 'elec_mobility_data.dat' 
file_hole= 'hole_mobility_data.dat'

f(x) = a1*exp(-b1/x)+c1*exp(-d1/x) #elec oop
g(x) = a2*exp(-b2/x)+c2*exp(-d2/x) #elec ip
h(x) = a3*exp(-b3/x)+c3*exp(-d3/x) #hole oop
j(x) = a4*exp(-b4/x)+c4*exp(-d4/x) #hole ip

fit f(x) file_elec u 4:(1/$5) via a1,b1,c1,d1
fit g(x) file_elec u 4:(1/$6) via a2,b2,c2,d2
fit h(x) file_hole u 4:(1/$5) via a3,b3,c3,d3
fit j(x) file_hole u 4:(1/$6) via a4,b4,c4,d4

plot f(x), file_elec u 4:(1/$5) w lp
plot g(x), file_elec u 4:(1/$6) w lp
plot h(x), file_hole u 4:(1/$5) w lp
plot j(x), file_hole u 4:(1/$6) w lp

