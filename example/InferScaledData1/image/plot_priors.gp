

gaussian(x,mean, precision) = sqrt(precision/(2.0*pi)) * exp(-0.5*precision*(x-mean)**2)
gammad(x, shape, iscale) = (1.0/gamma(shape))*(iscale**shape)*(x**(shape-1.0))*exp(-iscale*x)
set yrange [0:0.02]
set ytics 0.005

set term postscript eps enhanced color
set output "priors.eps"


set multiplot
set size  1.0, 0.5
set origin 0, 0.5
set xrange [-60:60]
plot gaussian(x,0.0,0.001) t "gaussian prior"

set origin 0,0
set xrange [0:60]
plot  gammad(x,1.0,0.01) t "gamma prior"
unset multiplot 

set term pop
set term png
set output "priors.png"

set multiplot
set size  1.0, 0.5
set origin 0, 0.5
set xrange [-60:60]
plot gaussian(x,0.0,0.001) t "prior for the mean"

set origin 0,0
set xrange [0:60]
plot  gammad(x,1.0,0.01) t "prior for the precision"
unset multiplot 


set term pop
set output