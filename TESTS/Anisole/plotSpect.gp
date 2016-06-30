set term png
set output "spectra.png"
set xrange [0:1900]
set yrange [0:0.02]
plot "shift.spe" u ($1-39284.7120):2 w l title "FC"
