reset
set term x11 size 680,620
set size 1.0, 1.0
set origin 0.0, 0.0
set yrange[-10.0:20.0]
plot 'bands_abinit_elk.dat' u 1:2 w lp title "ABINIT", 'bands_abinit_elk.dat' u 1:3 w lp title "Elk", 'bands_abinit_elk.dat' u 1:(0.0) w l title "Fermi level"



