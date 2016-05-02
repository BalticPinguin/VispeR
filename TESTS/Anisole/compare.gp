set term png
set output "compare.png"
plot "shift.spe" u ($1+39300):2 w l title "new",\
      "smsc.spe" u ($1+39300):2 w l title "old"

