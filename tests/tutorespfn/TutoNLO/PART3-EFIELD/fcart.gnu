set terminal postscript eps enhanced 22
set output 'fcart.eps'
set xrange [-5:5]
#set yrange [-0.012:0.01]
set ylabel 'Force (a.u.)'
set xlabel 'Electric field (10^{-4} a.u.)'
set size 1,1
set ticscale 2

set border 15 linewidth 1.5
set xtics 2
set mxtics 2

ucvol = 3.0109659E+02

f1(x) = a1*x**2 + b1*x + c1
f2(x) = a2*x**2 + b2*x + c2
f3(x) = a3*x**2 + b3*x + c3

fit f1(x) 'fcart.lst' u 1:2 via a1,b1,c1
fit f2(x) 'fcart.lst' u 1:2 via a2,b2,c2
fit f3(x) 'fcart.lst' u 1:2 via a3,b3,c3

plot "fcart.lst" u ($1*10000):($2) notitle w p pt 6 ps 2, \
   f1(x/10000) notitle w lines -1

print "Zstar = ", b1, b2, b3
print "dchi/dtau = ",a1/ucvol, a2/ucvol, a3/ucvol

set terminal x11
set size 1,1
replot
