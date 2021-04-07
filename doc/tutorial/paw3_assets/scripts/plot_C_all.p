reset
set term x11 size 680,620
set size 1.0, 1.0
set origin 0.0, 0.0
unset key
set multiplot
set yrange[-30:30]
set size 0.33, 0.33
set origin 0.0, 0.66
set title "logderiv.0"
plot "logderiv.0" u 1:2 w l, "logderiv.0" u 1:3 w l
set size 0.33, 0.33
set origin 0.33, 0.66
set title "logderiv.1"
plot "logderiv.1" u 1:2 w l, "logderiv.1" u 1:3 w l
set size 0.33, 0.33
set origin 0.66, 0.66
set title "logderiv.2"
plot "logderiv.2" u 1:2 w l, "logderiv.2" u 1:3 w l
set size 0.33, 0.33
set origin 0.00, 0.33
set yrange[-5:5]
set title "wfn (ae,ps,proj) for l=0"
plot 'wfn1' u 1:2 w l lc 1, 'wfn1' u 1:3 w l lc 2, 'wfn1' u 1:4 w l lc 3
set size 0.33, 0.33
set origin 0.33, 0.33
set title "wfn (ae,ps,proj) for l=1"
plot 'wfn2' u 1:2 w l lc 1, 'wfn2' u 1:3 w l lc 2, 'wfn2' u 1:4 w l lc 3
set size 0.33, 0.33
set origin 0.00, 0.0
set xrange[0:10]
set yrange[-3:5]
set title "vloc"
plot 'vloc' u 1:2 w l lc 1
set size 0.33, 0.33
set origin 0.66, 0.0
set xrange[0:10]
set yrange[-0.5:1]
set title "tprod"
plot 'tprod.1' u ($1*$1*0.5):2 w l lc 1, 'tprod.2' u ($1*$1*0.5):2 w l lc 2
unset multiplot

