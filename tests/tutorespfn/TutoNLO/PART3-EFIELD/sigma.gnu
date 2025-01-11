set terminal postscript eps enhanced 22
set output 'sigma.eps'
set xrange [-5:5]
#set yrange [-0.012:0.01]
set ylabel 'sigma (10^{-5} a.u.)'
set xlabel 'Electric field (10^{-4} a.u.)'
set size 1,1
set ticscale 2

set border 15 linewidth 1.5
set xtics 2
set mxtics 2

Ha_eV = 27.2113834
Bohr_Ang = 0.5291772083
e_Cb = 1.602176462e-19
HaBohr3_GPa=Ha_eV/Bohr_Ang**3*e_Cb*1.0e+21


f(x) = a*x + b 
fit f(x) 'sigma.lst' u 1:2 via a,b
plot "sigma.lst" u ($1*10000):($2*1e5) notitle w p pt 6 ps 2, \
   f(x/10000)*1e5 notitle w lines -1

print "piezoelectric constant = ", -a*e_Cb/(Bohr_Ang*1.0e-10)**2


set terminal x11
set size 1,1
replot
