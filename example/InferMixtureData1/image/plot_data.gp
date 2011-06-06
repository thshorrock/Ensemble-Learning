set size square
set size 0.6
set ytics 100
set xtics 1
set mxtics 2
set xrange [0:9]

set xlabel "data value"
set ylabel "sample"

set term postscript eps enhanced color
set output "data_points.eps"
plot "./data.dat" u 1:0 notitle
set term pop
set term png
set output "data_points.png"
replot
set term pop
set output