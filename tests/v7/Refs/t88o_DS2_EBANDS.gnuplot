# File to plot phonon bandstructure with gnuplot
#set terminal postscript eps enhanced color font 'Times-Roman,26' lw 2
#use the next lines to make a nice figure for a paper
#set term postscript enhanced eps color lw 0.5 dl 0.5
#set pointsize 0.275
set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )
unset key
# can make pointsize smaller (~0.5). Too small and nothing is printed
set pointsize 0.8
set view 0,0
set xrange [0:28]
set yrange [ -1.15248980E+01:  1.66482380E+01]
set xlabel "Momentum"
set ylabel "Energy [eV]"
set title "t88o\\_DS2\\_EBANDS.data"
set grid xtics
mband = 5
plot for [i=2:mband] "t88o_DS2_EBANDS.data" u 1:i every :1 with lines linetype -1
pause -1
