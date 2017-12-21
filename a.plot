#!/usr/bin/gnuplot
# this is a comment
reset
# the next line tells gnuplot to keep the plot window open
# the rest is the same as before
set term postscript enhanced
set output "py_c.eps"
set title "Intel(R) Xeon(R) CPU E5-2640 v4 @ 2.40GHz"
set xlabel "   "
set ylabel "Time(ms)"
set logscale y
set xrange [-0.5:2.5]
set style line 1 lt 1 lc 1 lw 2 pt 1 pi -1 ps 1.5
set style line 2 lt 1 lc 3 lw 2 pt 1 pi -1 ps 1.5
set style line 3 lt 4 lc 1 lw 2 pt 1 pi -1 ps 1.5
set style line 4 lt 4 lc 3 lw 2 pt 1 pi -1 ps 1.5

set style fill solid
set boxwidth 0.5
set xtics ("Python" 0.25, "C++" 1.75)

plot "a.dat" every 2 using 1:2 with boxes ls 1 title "DFS", \
     "a.dat" every 2::1 using 1:2 with boxes ls 2 title "Neighboring List"

