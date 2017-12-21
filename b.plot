#!/usr/bin/gnuplot
# this is a comment
reset
# the next line tells gnuplot to keep the plot window open
# the rest is the same as before
set term postscript enhanced
set output "openmp.eps"
set title "Intel(R) Xeon(R) CPU E5-2640 v4 @ 2.40GHz"
set yrange [0:60]
set xlabel "N"
set ylabel "Time(ms)"
set style line 1 lt 1 lc 1 lw 2 pt 1 pi -1 ps 1.5
set style line 2 lt 1 lc 3 lw 2 pt 1 pi -1 ps 1.5
set style line 3 lt 4 lc 1 lw 2 pt 1 pi -1 ps 1.5
set style line 4 lt 4 lc 3 lw 2 pt 1 pi -1 ps 1.5


plot "b.dat" u 1:2 with linespoints ls 1 title "DFS", \
     "b.dat" u 1:3 with linespoints ls 2 title "Neighboring List"

