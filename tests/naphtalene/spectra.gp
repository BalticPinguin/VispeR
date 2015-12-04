set yrange [0:3]

plot "DR_full.spe" u 1:2 w l title "DR", "shift.spe" w l title "FC", "changed.spe" u 1:2 w l title "changed"

pause -1
