set xrange [17500:21500]
set key left

set term png
set output "models.png"

plot "shift_stretch.spe" u 1:2 w l title "shift",\
"changed.spe" w l title "changed",\
"Duschinsky.spe" w l title "Duschinsky, OPA",\
"Duschinsky_full.spe" w l title "Duschinsky, full"
