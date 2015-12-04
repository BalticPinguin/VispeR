set xrange [00:1800]
set yrange [-0.01:0.4]
#plot "DR.spe" w l title "DR", "smsc.spe" w l title "smsc", "changed2.spe" w l title "changed", "shift2.spe" u 1:($2) w l title "shift", "DR.spe" w l title "DR2"
#pause -1

#unset key
#set pm3d map
#set pm3d interpolate 2,2
#set xrange [0:42]
#set yrange [0:42]
#splot "J.mat" matrix
#pause -1

set term png
set output "my_spects.png"

set notics
#set xtics 200 ""
set multiplot

set origin 0.00,0.70
set size 1.00, 0.33
plot "DR.spe" u ($1+39300):2 w l title "Duschinsky"

set origin 0.00,0.45
set size 1.00, 0.33
plot "smsc.spe" u ($1+39300):2 w l title "broadened"

set origin 0.00,0.20
set size 1.00, 0.33
plot "changed.spe" u ($1+39300):2 w l title "changed"

set origin 0.00,-0.04
set size 1.00, 0.33
plot "shift.spe" u ($1+39300):2 w l title "shift"

#plot "G09-spectrum.spe" u ($1+39300):2 w l
#pause -1

unset multiplot

unset xrange
set output "compare.png"
plot "shift_11_15_old.spe" u ($1+39300):2 w l title "old",\
      "shift_11_15.spe" u ($1+39300):2 w l title "new"
