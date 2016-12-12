# File to plot bandstructure with gnuplot
set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )
unset ztics
unset key
# can make pointsize smaller (~0.5). Too small and nothing is printed
set pointsize 0.8
set grid xtics
set view 0,0
set xrange [0:40]
set yrange [  0.00000000E+00:  3.89492578E+01]
#use the next lines to make a nice figure for a paper
#set term postscript enhanced eps color lw 0.5 dl 0.5
#set pointsize 0.275
nbranch = 3
plot for [i=2:nbranch] "t88o_DS2_PHBANDS.data" u 1:i every :1 with lines linetype -1
pause -1
