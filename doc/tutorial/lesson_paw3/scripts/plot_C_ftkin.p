reset
set xrange[0:20]
plot 'ftkin.1' u ($1*$1*0.5):2 w l lc 1, 'ftkin.2' u ($1*$1*0.5):2 w l lc 2

