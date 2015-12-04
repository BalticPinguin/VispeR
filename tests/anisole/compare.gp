set term png
set output "compare.png"
plot "shift_11_15_old.spe" u ($1+39300):2 w l title "old",\
      "shift_11_15.spe" u ($1+39300):2 w l title "new"

