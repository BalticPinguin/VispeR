set xrange [-43000:-36000]
set yrange [-0.01:1.1]
plot "changed2.spe" w l title "changed", "changed2R.spe" u ($1-2*39284.7120575):2  w l title "changed_R", "changed2R.log" u ($2-2*39284.7120575):1  w p title " "
pause -1
