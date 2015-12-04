plot "spectrum.dat" u 1:2 w p title "reference, LIF",\
     "spectrum.dat" u 1:($3*10) w p title "reference, ICLAS"
#   "DR.spe" w l title "DR",\
#   "changed.spe" w l title "changed",\
#   "shift.spe" w l title "shift",\
#   "gradient.spe" w l title "gradient"
pause -1
