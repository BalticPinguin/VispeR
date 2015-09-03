#set xrange [22000:34000]
set yrange [0:3]
#plot "DR.spe" u 1:2 w l title "DR_changed", "shift.spe" w l title "FC", "changed.spe" u 1:2 w l title "changed", "DR_same.spe" u 1:2 w l title "DR_same"
#plot "shift.spe" w l title "FC", "DR_same.spe" u 1:2 w l title "DR_same", "DR" u ($5+32000):4:1 w p palette title " ","shift" u ($5+32000):4:1 w p palette title " "

plot "DR.spe" u 1:2 w l title "DR_changed", "shift.spe" w l title "FC", "changed.spe" u 1:2 w l title "changed", "DR_same.spe" u 1:2 w l title "DR_same", "DR_full.spe" u 1:2 w l title "DR_full"

#plot "DR.spe" u 1:2 w l title "DR_changed", "changed.spe" u 1:2 w l title "changed", "DR_full.spe" u 1:2 w l title "DR_full", "DR" u 2:1 w p title " "
pause -1
