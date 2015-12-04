#set xrange [22000:34000]
set yrange [0:5]
#plot "DR.spe" u 1:2 w l title "DR_changed", "shift.spe" w l title "FC", "changed.spe" u 1:2 w l title "changed", "DR_same.spe" u 1:2 w l title "DR_same"
#plot "shift.spe" w l title "FC", "DR_same.spe" u 1:2 w l title "DR_same", "DR" u ($5+32000):4:1 w p palette title " ","shift" u ($5+32000):4:1 w p palette title " "
#plot "DR.spe" u 1:2 w l title "DR_changed", "changed.spe" u 1:2 w l title "changed", "DR_full.spe" u 1:2 w l title "DR_full", "DR" u 2:1 w p title " "

set term png
set output "spectra.png"
set xrange [22000:34000]

set multiplot
set origin 0.,0.
set size 1.0, 1.0
plot "DR_changed.spe" u 1:2 w i title "DR_changed", "shift.spe" w i title "FC", "changed.spe" u 1:2 w i title "changed", "DR_same.spe" u 1:2 w l title "DR_same", "DR_full.spe" u 1:2 w l title "DR_full"

set origin 0.02, 0.6
set size 0.5, 0.4
set xrange [30000:33000]
set yrange [0:0.6]
unset xtics
unset ytics
unset key
plot "DR_changed.spe" u 1:2 w i title "DR_changed", "shift.spe" w i title "FC", "changed.spe" u 1:2 w i title "changed", "DR_same.spe" u 1:2 w l title "DR_same", "DR_full.spe" u 1:2 w l title "DR_full"

pause -1
