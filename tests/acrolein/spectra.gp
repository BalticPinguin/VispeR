set term png
set output "compare2.png"
set xrange [17500:21500]
set key left
#plot "t1-s0-broad.spe" u ($1-146):($2/280) w l title "g09",\
#"tests/shift.spe" u 1:($2*0.95) w l title "shift",\
#"tests/fit.spe" u 1:($2*0.95) w l title "fitted"

#set xlabel "HR fit"
#set ylabel "HR"
#set logscale x
#set logscale y
#set yrange [0.02:1.2]
#set xrange [0.02:1]
#plot "tests/compare.dat" u 1:3:($2-$4) w p pt 7 palette notitle, 4/pi*x notitle

plot "t1-s0-broad.spe" u ($1-146):($2/280) w l title "g09",\
"tests/fit.spe" u 1:($2*0.95) w l title "fitted",\
"tests/shift_stretch.spe" u 1:2 w l title "shift scaled",\

#plot "tests/shift.spe" u 1:($2*0.95) w l title "shift",\
#"Tobias/shift.spe" u 1:($2*0.95) w l title "shift",\
#"Tobias/Duschinsky_full.spe" u 1:2 w l title "DR"
#pause -1



set output "compare3.png"
plot "tests/shift_stretch.spe" u 1:2 w l title "shift scaled",\
     "tests/shift_11_15.spe" u 1:2 w l title "shift new"
