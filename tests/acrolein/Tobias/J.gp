set term png
set output "J.png"

unset key

set pm3d map
set pm3d interpolate 2,2
set xrange [0:17]
set yrange [0:17]

set multiplot

set origin 0.1, 0.0
set size 0.9, 1.0
splot "J.mat" matrix

unset xtics
unset ytics
unset cbtics
#unset cb
set xrange [0:1]
set origin 0.0, 0.0
set size 0.2, 1.0
splot "K.mat" matrix
