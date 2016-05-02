set term png

#set output "FC_fluor.png"
#set xrange [34500:42000]
#plot "FC_fluor.spe" u ($1+40000):($2) w l t "FC", "experiment.dat" u (10000000/$1):($2/10) w l t "experiment"

set output "FC_abs.png"
#set xrange [34500:34800]
set xrange[0:3000]
plot "FC_abs.spe" u ($1-1256):($2) w l t "FC" , "experiment.dat" u (10000000/$1-36150):($2/3) w l t "experiment"

set output "the_abs.png"
plot "DR_abs.spe" u ($1-1256):($2*8) w l t "DR", "FC_abs.spe" u ($1-1256):($2) w l t "FC" , "experiment.dat" u (10000000/$1-36150):($2/3) w l t "experiment"
