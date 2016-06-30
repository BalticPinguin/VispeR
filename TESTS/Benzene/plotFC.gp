set term png
set output "FC_fluor.png"
plot "FC_fluor.spe" u ($1):($2*10000) w l t "FC", "experiment_fluer.dat" u (10000000/$1-35000):($2/16000) w l t "experiment"

set output "FC_abs.png"
#set xrange [41000:47000]
plot "FC_exc.spe" u ($1):2 w l t "FC", "exp1.dat" u ($1-35000):($2/60) w l t "experiment1"

set output "FC_abs2.png"
unset xrange
#set xrange [40000:47000]
plot "FC_exc.spe" u ($1):2 w l t "FC (singlet)"

set output "abs_FCvsDR.png"
plot "FC_exc.spe" u ($1):2 w l t "FC","DR_exc.spe" u ($1):($2*1) w l t "DR"
