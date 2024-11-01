# File to plot electron bandstructure with gnuplot
#use the next lines to make a nice figure for a paper
#set term postscript enhanced eps color lw 0.5 dl 0.5
#set pointsize 0.275
set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )
unset key
# can make pointsize smaller (~0.5). Too small and nothing is printed
set pointsize 0.8
set view 0,0
set xrange [0:40]
set yrange [ -4.06419318E+00:  3.89483730E+01]
set xlabel "Momentum"
set ylabel "Energy [meV]"
set title "t88o\\_DS2\\_PHBANDS.data"
# Add vertical lines in correspondence of high-symmetry points.
unset xtics
set arrow from 0,graph(0,0) to 0,graph(1,1) nohead
set arrow from 20,graph(0,0) to 20,graph(1,1) nohead
set arrow from 40,graph(0,0) to 40,graph(1,1) nohead
nbranch = 3
plot for [i=2:nbranch] "t88o_DS2_PHBANDS.data" u 1:i every :1 with lines linetype -1
pause -1
