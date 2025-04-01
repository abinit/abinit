set terminal postscript eps enhanced 22
set output 'pol.eps'
set xrange [-5:5]
#set yrange [-0.012:0.01]
set ylabel 'Polarization (a.u.)'
set xlabel 'Electric field (10^{-4} a.u.)'
set size 1,1
set ticscale 2

set border 15 linewidth 1.5
set xtics 2
set mxtics 2
set key 0,-0.0262

eps0     = 8.854187817e-12
e_Cb     = 1.602176462e-19
Bohr_Ang = 0.5291772083

f1(x) = 2*a1*x**2 + b1*x + c1
f2(x) = 2*a2*x**2 + b2*x + c2
f3(x) = 2*a3*x**2 + b3*x + c3
fit f1(x) 'pol.lst' u 1:2 via a1,b1,c1
fit f2(x) 'pol.lst' u 1:3 via a2,b2,c2
fit f3(x) 'pol.lst' u 1:4 via a3,b3,c3

plot "pol.lst" u ($1*10000):($2) title "P_x" w p pt 6 ps 2, \
     "pol.lst" u ($1*10000):($3) title "P_y" w p pt 8 ps 2, \
     "pol.lst" u ($1*10000):($4) title "P_z" w p pt 4 ps 2, \
   f1(x*0.0001) notitle w lines -1, \
   f2(x*0.0001) notitle w lines -1, \
   f3(x*0.0001) notitle w lines -1

print "Dielectric constant (P_x) = ", (1 + 4*pi*b1)
print "Dielectric constant (P_y) = ", (1 + 4*pi*b2)

chi2 = 4*pi*a3*(Bohr_Ang**2)*1.0e-8*4*pi*eps0/e_Cb
print "Non-linear optical susceptibility (pm/V) d = ", chi2/2

set terminal x11
set size 1,1
replot
