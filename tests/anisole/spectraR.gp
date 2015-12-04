set xrange [-43000:-36000]
set yrange [-0.01:1.1]
plot "DR.spe" w l title "DR", "DRR.spe" u ($1-2*39284.7120575):2 w l title "DR_R",\
   "DR_opa.spe" w l title "DR opa", "DRR_opa.spe" u ($1-2*39284.7120575):2 w l title "DR_R opa"
    #"smsc.spe" w l title "smsc", "smscR.spe" u ($1-2*39284.7120575):2 w l title "smsc_R",\
    #"changed2.spe" w l title "changed", "changed2R.spe" u ($1-2*39284.7120575):2  w l title "changed_R",\
    #"shift2.spe" u 1:($2) w l title "shift", "shift2R.spe" u ($1-2*39284.7120575):2 w l title "shift_R"
pause -1
